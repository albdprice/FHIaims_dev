!****s* FHI-aims/integrate_hessian
!  NAME
!   integrate_hessian
!  SYNOPSIS

subroutine integrate_hessian &
     ( hartree_potential_std,v_hartree_gradient_std, & 
       first_order_potential_moving_grid_std, & 
       rho_std,rho_gradient_std,  &
       partition_tab_std, partition_deriv_std, basis_l_max,     &
       KS_ev_compute, & 
       density_matrix, first_order_density_matrix,  & 
       energy_density_matrix, first_order_energy_density_matrix,  &  
       max_occ_number, &
       hessian,gradient_dipole & 
     ) 

!  PURPOSE
!  Integrates the hessian(alpha,beta)
!  (0) changed from integrate_first_order_H 
!  (1) and then add hessian_wave part
!  (2) and calc hessian_hellman_feynman,H_pulay1,H_pulay2, S1,S2 at every bathes
!  (3) after all bathes, get hessian_pulay. 
!  
!  shanghui 2013.02.21 @ Berlin 
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data 
  use load_balancing
  use cartesian_ylm  !shanghui add for Hessian
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: hartree_potential_std
  real*8, target, dimension(3, n_atoms, n_full_points) ::v_hartree_gradient_std
  real*8, target, dimension(3, n_atoms, n_full_points) :: & 
                            first_order_potential_moving_grid_std
  real*8, target, dimension(n_spin, n_full_points)    :: rho_std
  real*8, target, dimension(3, n_spin, n_full_points) :: rho_gradient_std
  real*8, target, dimension(n_full_points)            :: partition_tab_std
  real*8, target, dimension(3,n_atoms,n_full_points)  :: partition_deriv_std

  integer ::  basis_l_max (n_species)
  real*8,   dimension(n_basis, n_states), intent(in) :: KS_ev_compute

  real*8, dimension(n_basis,n_basis), intent(in) :: density_matrix
  real*8, dimension(n_basis,n_basis), intent(in) :: energy_density_matrix
  real*8, dimension(3, n_atoms, n_basis,n_basis),intent(in):: first_order_density_matrix
  real*8, dimension(3, n_atoms, n_basis,n_basis),intent(in):: first_order_energy_density_matrix
  integer, dimension(n_spin), intent(in) :: max_occ_number
 
  real*8, dimension(3*n_atoms, 3*n_atoms), intent(out) :: hessian
  real*8, dimension(3*n_atoms, 3), intent(out) ::  gradient_dipole


!  INPUTS
!  o hartree_potential_std -- Hartree potential
!  o rho_std -- electron density
!  o rho_gradient_std -- gradient of electron density.
!    These should only ever be referenced if (use_gga)
!    import dimensions from above (if not used, all dimensions=1)
!  o partition_tab_std -- values of partition functions
!  o basis_l_max -- maximum l of basis functions.
!
!  OUTPUT
!  o hellman_feynman_hessian 
!  o pulay_hessian
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



  !  local variables

  real*8, dimension(3, n_atoms, 3,n_atoms) :: hellman_feynman_hessian
  real*8, dimension(3, n_atoms, 3,n_atoms) :: pulay_hessian
  real*8, dimension(n_spin) :: local_potential_parts

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

  real*8 coord_current(3)

!  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
!  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)

  real*8,dimension(:,:),allocatable:: dist_tab
  real*8,dimension(:,:),allocatable:: dist_tab_sq

  real*8 i_r(n_max_compute_atoms)

!  real*8 dir_tab(3,n_centers_integrals, n_max_batch_size)
  real*8, dimension(:,:,:),allocatable:: dir_tab
  real*8, dimension(:,:,:),allocatable:: dir_tab_backup
  real*8,dimension(:,:),allocatable:: dist_tab_backup
  real*8,dimension(:,:),allocatable:: dist_tab_sq_backup
 

  real*8 trigonom_tab(4,n_max_compute_atoms)

  real*8,dimension(:,:,:),allocatable:: H_times_psi
  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: kinetic_wave
  real*8,dimension(:,:)  ,allocatable:: wave
  real*8, dimension( :, :, :), allocatable :: local_first_order_rho
  real*8, dimension( :, :, :), allocatable :: local_first_order_rho_moving_grid
  real*8, dimension( :, :, :), allocatable :: first_order_DM_rho


  real*8, dimension(:),    allocatable :: en_density_xc
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8, dimension(:, :), allocatable :: local_xc_derivs
  real*8, dimension(:, :), allocatable :: local_dVxc_drho
  real*8, dimension(:,:,:),allocatable :: xc_gradient_deriv

  real*8, dimension(:,:),  allocatable :: local_rho
  real*8, dimension(:,:,:),allocatable :: local_rho_gradient


  real*8, dimension(:,:,:),  allocatable :: local_v_hartree_gradient
  real*8, dimension(:,:,:),  allocatable :: local_first_order_potential_moving_grid

  !     Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
  !     The hope is that such a separate treatment will allow to minimize numerical noise
  !     introduced through ZORA
  real*8, dimension(:), allocatable :: hamiltonian_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)
  real*8 :: partition_deriv_local(3,n_atoms,n_max_batch_size) 
  real*8 :: energy_partition(n_max_batch_size)

  real*8, dimension(:,:), allocatable :: gradient_basis_wave
  real*8, dimension(:,:,:), allocatable :: gradient_basis_wave_npoints

!-------For pulay_hessian--------------------------------------
    real*8, dimension(:,:,:), allocatable :: hessian_basis_wave
    real*8, dimension(:,:,:), allocatable :: sum_hessian
    real*8, dimension(:,:,:), allocatable :: cartesians
    real*8, dimension(:,:,:), allocatable :: sum_gradient
    real*8,dimension(:),allocatable:: radial_wave_2nd_deriv
    real*8, dimension(:,:,:,:), allocatable :: H_times_gradient_psi
    real*8,dimension(:)  ,allocatable:: kinetic_wave_deriv
    real*8, dimension(:,:,:), allocatable :: kinetic_gradient_basis_wave

    real*8, dimension(:,:,:,:), allocatable :: first_order_S
    real*8, dimension(:,:,:,:), allocatable :: first_order_H_pulay
    real*8, dimension(:,:,:,:,:,:), allocatable :: second_order_S
    real*8, dimension(:,:,:,:,:,:), allocatable :: second_order_H_pulay
         
    real*8  pulay_force(3,n_atoms)
!-------End for pulay_hessian----------------------------------

!-------For ionic_hessian-------------------------------------
   real*8   ionic_hessian(3,n_atoms,3,n_atoms)
   real*8   direction(3), distance_squared 
!-------End for ionic_hessian-------------------------------------

!-------For all HF and Pulay Hessian's moving point---------------
   integer :: atom_point(n_max_batch_size)
!-------End for all HF and Pulay Hessian's moving point---------------
   integer :: current_atom



  !     Following is all that is needed for the handling of ZORA scalar relativity

  real*8, dimension(n_spin) :: zora_operator
  logical, dimension(n_spin) :: t_zora
  real*8, dimension(n_spin) :: zora_potential_parts

  real*8, dimension(:), allocatable :: dist_tab_full
  real*8, dimension(:,:), allocatable :: dir_tab_full_norm
  real*8, dimension(:), allocatable :: i_r_full

  real*8, dimension(:,:,:,:), allocatable :: zora_vector1
  real*8, dimension(:,:,:,:), allocatable :: zora_vector2

  ! This term contains contributions from the xc potential and the
  ! zora formalism (if applicable) which are summed up using Gauss' law:
  ! < grad(phi_i) | local_gradient_sum |grad(phi_j) >
  real*8, dimension(3,n_spin) :: sum_of_local_gradients

  !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

  integer :: n_compute_c, n_compute_a
!  integer :: i_basis(n_centers_basis_I)
  integer,dimension(:),allocatable :: i_basis

  integer :: n_compute_fns

!  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
!  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
!  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer,dimension(:),  allocatable :: i_basis_fns
  integer,dimension(:,:),allocatable :: i_basis_fns_inv
  integer,dimension(:),  allocatable :: i_atom_fns

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms)

  ! indices for basis functions that are nonzero at current point

  integer :: rad_index(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_ham)
  integer :: l_index(n_max_compute_fns_ham)
  integer :: l_count(n_max_compute_fns_ham)
  integer :: fn_atom(n_max_compute_fns_ham)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute
  integer :: zero_index_point(n_max_compute_ham)

  ! active atoms in current batch
  integer :: n_batch_centers
  integer :: batch_center(n_centers_integrals)

  !     for splitting of angular shells into "octants"

  integer division_low
  integer division_high

  !  counters

  integer i_basis_1, j_basis
  integer i_atom, i_atom_2 ,j_atom
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord, j_coord,i_coord_1, i_coord_2
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2
  integer :: i_full_points_first_order_rho  
 
  integer i_compute 
 
  integer i_counter
  integer index_hessian(3,3)

  integer :: i_spin
  character*200 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used


  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: partition_deriv(:,:,:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: rho_gradient(:,:,:)
  real*8, pointer :: hartree_potential(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i_off, i, j, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

   !----------shanghui add here for evaluate_first_order_rho----
   integer n_occ_states
    n_occ_states= max_occ_number(1)  ! now only for nspin=1 case, in the futher, need extended
   !----------shanghui end add here for evaluate_first_order_rho----

  ! begin work


  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating Hessian matrix: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating Hessian matrix: batch-based integration."
  endif
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

  ! begin with general allocations

  allocate(dist_tab(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab                      ')
  allocate(dist_tab_backup(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab_backup                      ')

  allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab_sq                   ')
  allocate(dist_tab_sq_backup(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab_sq_backup                      ')

  allocate(dir_tab(3,n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dir_tab                       ')

  allocate(dir_tab_backup(3,n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dir_tab_backup                       ')


  allocate(i_basis_fns(n_basis_fns*n_centers_integrals), stat=info)
  call check_allocation(info, 'i_basis_fns                   ')

  allocate(i_basis_fns_inv(n_basis_fns,n_centers), stat=info)
  call check_allocation(info, 'i_basis_fns_inv               ')

  allocate(i_atom_fns(n_basis_fns*n_centers_integrals),stat=info)
  call check_allocation(info, 'i_atom_fns                    ')

  allocate( en_density_xc(n_max_batch_size),stat=info)
  call check_allocation(info, 'en_density_xc                 ')

  allocate( local_xc_derivs(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_xc_derivs               ')
  allocate( local_dVxc_drho(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_dVxc_drho            ')

  allocate( xc_gradient_deriv(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'xc_gradient_deriv             ')

  allocate( local_rho(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho                     ')


  allocate( local_rho_gradient(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho_gradient            ')

  allocate( local_v_hartree_gradient(3,n_atoms,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_v_hartree_gradient                     ')

  allocate( local_first_order_potential_moving_grid(3,n_atoms,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_first_order_potential_moving_grid     ')

!shanghui-gga: we do not gga, but we always need gradient. 
!  if ((flag_rel.eq.REL_none.or.flag_rel==REL_atomic_zora.or.flag_rel.eq.REL_own).and.(.not.(use_gga))) then
     !       no gradients needed
!     l_ylm_max = l_wave_max
!  else if ((flag_rel.eq.REL_zora).or.(use_gga).or.(flag_rel==REL_KOLNING_HARMON)) then
     l_ylm_max = l_wave_max
     allocate (gradient_basis_wave(n_max_compute_ham,3),STAT=info)
     call check_allocation(info, 'gradient_basis_wave           ')
     allocate (gradient_basis_wave_npoints(n_max_compute_ham,3,n_max_batch_size),STAT=info)
     call check_allocation(info, 'gradient_basis_wave_npoints        ')

     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')

!  end if !-----shanghui always calcualet gradients-------------------------

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info )
  call check_allocation(info, 'index_lm                      ')

  allocate ( hamiltonian_shell(n_max_compute_ham*n_max_compute_ham), STAT=info )
  call check_allocation(info, 'hamiltonian_shell             ')

  allocate(H_times_psi(n_max_compute_ham, n_max_batch_size, n_spin), STAT=info )
  call check_allocation(info, 'H_times_psi                   ')

  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(kinetic_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'kinetic_wave                  ')

  allocate(wave(n_max_compute_ham, n_max_batch_size), STAT=info )
  call check_allocation(info, 'wave                          ')

  allocate(local_first_order_rho(3,n_atoms, n_max_batch_size), STAT=info )
  call check_allocation(info, 'local_first_order_rho                          ')
  allocate(local_first_order_rho_moving_grid(3,n_atoms, n_max_batch_size), STAT=info )
  call check_allocation(info, 'local_first_order_rho_moving_grid               ')

  allocate(first_order_DM_rho(3,n_atoms, n_max_batch_size), STAT=info )
  call check_allocation(info, 'first_order_DM_rho                          ')

  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')


!-------For pulay_hessian--------------------------------------
    if (.not.allocated(cartesians)) then
       allocate(cartesians(n_max_cartesian, 0:l_wave_max,n_max_compute_atoms ),stat=info)
       call check_allocation(info, 'cartesians                    ')
    end if
    if (.not.allocated(sum_gradient)) then
       allocate(sum_gradient(3, (l_wave_max+1) ** 2, n_max_compute_atoms ),stat=info)
       call check_allocation(info, 'sum_gradient                  ')
    end if
    if (.not.allocated(sum_hessian)) then
       allocate(sum_hessian(6, (l_wave_max+1) ** 2, n_max_compute_atoms ),stat=info)
       call check_allocation(info, 'sum_hessian                   ')
    end if
    if (.not.allocated(hessian_basis_wave)) then
       allocate(hessian_basis_wave(n_max_compute_ham, 6, n_max_batch_size),stat=info)
       call check_allocation(info, 'hessian_basis_wave            ')
    end if
    if(.not. allocated(radial_wave_2nd_deriv))then
       allocate(radial_wave_2nd_deriv(n_max_compute_fns_ham),stat=info)
       call check_allocation(info, 'radial_wave_2nd_deriv         ')
    end if
    if(.not.allocated(H_times_gradient_psi)) then
       allocate(H_times_gradient_psi(n_max_compute_ham,3, n_max_batch_size, n_spin), STAT=info )
       call check_allocation(info, 'H_times_gradient_psi                   ')
    endif
    allocate (kinetic_gradient_basis_wave(n_max_compute_ham,3,n_max_batch_size ),stat=info)
    call check_allocation(info, 'kinetic_gradient_basis_wave   ')

    allocate(kinetic_wave_deriv(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'kinetic_wave_deriv            ')

    allocate(first_order_S(3,n_atoms,n_basis,n_basis),stat=info)
    call check_allocation(info, 'first_order_S            ')
    allocate(second_order_S(3,n_atoms,3,n_atoms,n_basis,n_basis),stat=info)
    call check_allocation(info, 'second_order_S            ')

    !allocate(first_order_H(3,n_atoms,n_basis,n_basis),stat=info)
    !call check_allocation(info, 'first_order_H           ')

    allocate(first_order_H_pulay(3,n_atoms,n_basis,n_basis),stat=info)
    call check_allocation(info, 'first_order_H_pulay            ')
    allocate(second_order_H_pulay(3,n_atoms,3,n_atoms,n_basis,n_basis),stat=info)
    call check_allocation(info, 'second_order_H_pulay            ')
!-------End for pulay_hessian--------------------------------------




  if (flag_rel.eq.REL_zora.or.flag_rel==REL_KOLNING_HARMON ) then
     ! allocate all arrays relevant for ZORA

     if (.not.allocated(dist_tab_full)) then
        allocate(dist_tab_full(n_centers_integrals),STAT=info )
        call check_allocation(info, 'dist_tab_full                 ')

     end if
     if (.not.allocated(dir_tab_full_norm)) then
        allocate(dir_tab_full_norm(3,n_centers_integrals),STAT=info )
        call check_allocation(info, 'dir_tab_full_norm             ')
     end if
     if (.not.allocated(i_r_full)) then
        allocate(i_r_full(n_centers_integrals),STAT=info )
        call check_allocation(info, 'i_r_full                      ')

     end if

     if (.not.allocated(zora_vector1)) then
        allocate(zora_vector1(n_max_compute_ham,3,n_max_batch_size,n_spin),STAT=info )
        call check_allocation(info, 'zora_vector1                  ')
     end if
     if (.not.allocated(zora_vector2)) then
        allocate(zora_vector2(n_max_compute_ham,3,n_max_batch_size,n_spin),STAT=info )
        call check_allocation(info, 'zora_vector2                  ')
     end if

  end if

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(rho(n_spin,batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    allocate(hartree_potential(batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,1,hartree_potential_std,hartree_potential)

    if(use_density_gradient) then
      allocate(rho_gradient(3,n_spin,batch_perm(n_bp)%n_full_points))
      call permute_point_array(n_bp,3*n_spin,rho_gradient_std,rho_gradient)
    else
      ! Even though rho_gradient_std is allocated to a dummy size in this case,
      ! the array rho_gradient is used below as a dummy argument in full size
      ! (calls to evaluate_xc).
      ! rho_gradient therefore shouldn't be a dangling or nullified pointer
      ! since this will generated errors when in bounds checking mode.
      ! So we have to allocate it here, although it isn't needed actually.
      allocate(rho_gradient(3,n_spin,batch_perm(n_bp)%n_full_points))
    endif

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))


  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    partition_deriv => partition_deriv_std
    rho => rho_std
    hartree_potential => hartree_potential_std
    rho_gradient => rho_gradient_std

  endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------

  ! initialize


  i_basis_fns_inv = 0


  ! initialize index_lm

  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo


  !--------shanghui add for hessian_index-----------------
    i_counter = 0

    ! Since not all compilers initialize variables to zero index_hessian needs
    ! to be initialized here. Otherwise the max-construct below will potentially
    ! fill the array with nonsense. This happens at least with gfortran.
    index_hessian = 0

    do i_coord_1 = 1,3
       i_counter = i_counter +1
       index_hessian(i_coord_1, i_coord_1) = i_counter
       do i_coord_2 = i_coord_1+1,3
          i_counter = i_counter +1
          index_hessian(i_coord_1, i_coord_2) = i_counter
       end do
    end do

    do i_coord_1 = 1,3
       do i_coord_2 = 1,3
          index_hessian(i_coord_1, i_coord_2) = max(index_hessian(i_coord_1, i_coord_2), index_hessian(i_coord_2, i_coord_1))
       end do
    end do
  !--------shanghui end add for hessian_index--------------


  !--------shanghui add for ionic_hessian------------------------

     ionic_hessian(1:3,1:n_atoms, 1:3, 1:n_atoms) = 0.0d0


    ! if ( n_periodic == 0 .and. (.not. force_new_functional) ) then
       do i_atom = 1, n_atoms
       do j_atom = 1, n_atoms
           if (i_atom .ne. j_atom) then
                 distance_squared = 0.d0
                 do i_coord = 1, 3, 1
                    direction(i_coord) =   &
                         coords(i_coord, j_atom)-coords(i_coord, i_atom)
                    distance_squared = distance_squared +   &
                         direction(i_coord) * direction(i_coord)
                 end do

                 do i_coord = 1, 3, 1 
                 do j_coord = 1, 3, 1 
                    if(j_coord.eq. i_coord) then 

                    ionic_hessian(i_coord, i_atom, i_coord, j_atom) =   &
                    species_z(species(i_atom)) * species_z(species(j_atom))   &
                    /(distance_squared ** 1.5d0)                        & 
                    -3.0d0*direction(i_coord)*direction(j_coord)        & 
                    *species_z(species(i_atom)) * species_z(species(j_atom))   &
                    /(distance_squared ** 2.5d0)

                   !write(use_unit,*) i_coord, i_atom, i_coord, j_atom , & 
                   !ionic_hessian(i_coord, i_atom, i_coord, j_atom)
             
                    ionic_hessian(i_coord, i_atom, i_coord, i_atom) =   &
                    ionic_hessian(i_coord, i_atom, i_coord, i_atom) -   &
                    species_z(species(i_atom)) * species_z(species(j_atom))   &
                    /(distance_squared ** 1.5d0)                        & 
                    +3.0d0*direction(i_coord)*direction(j_coord)        & 
                    *species_z(species(i_atom)) * species_z(species(j_atom))   &
                    /(distance_squared ** 2.5d0)

                    else 

                    ionic_hessian(i_coord, i_atom, j_coord, j_atom) =   &
                    -3.0d0*direction(i_coord)*direction(j_coord)        & 
                    *species_z(species(i_atom)) * species_z(species(j_atom))   &
                    /(distance_squared ** 2.5d0)

                    ionic_hessian(i_coord, i_atom, j_coord, i_atom) =   &
                    ionic_hessian(i_coord, i_atom, j_coord, i_atom)    &
                    +3.0d0*direction(i_coord)*direction(j_coord)        & 
                    *species_z(species(i_atom)) * species_z(species(j_atom))   &
                    /(distance_squared ** 2.5d0)
                   
                    endif
                 enddo ! i_coord
                 enddo ! j_coord
        
           
            end if !i_atom .ne. j_atom
        enddo  ! i_atom
        enddo  ! j_atom
    
     !end if

  !--------shanghui end add for ionic_hessian------------------------

  hellman_feynman_hessian(1:3,1:n_atoms,1:3,1:n_atoms)=0.0d0
  

  first_order_H_pulay(1:3,1:n_atoms,1:n_basis,1:n_basis) = 0.0d0 
  first_order_S(1:3,1:n_atoms,1:n_basis,1:n_basis) = 0.0d0
  second_order_H_pulay(1:3,1:n_atoms,1:3,1:n_atoms,1:n_basis,1:n_basis) = 0.0d0
  second_order_S(1:3,1:n_atoms,1:3,1:n_atoms,1:n_basis,1:n_basis) = 0.0d0

  gradient_dipole(1:3*n_atoms,1:3) = 0.0d0

  i_full_points = 0
  i_full_points_2 = 0
  i_full_points_first_order_rho = 0


  ! perform partitioned integration, batch by batch of integration point.
  ! This will be the outermost loop, to save evaluations of the potential.
  ! and the Y_lm functions

  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()


  do i_my_batch = 1, n_my_batches_work, 1
     if(get_batch_weights) time_start = mpi_wtime()

     n_compute_c = 0
     n_compute_a = 0
     i_basis = 0

     i_point = 0

     ! loop over one batch
     do i_index = 1, batches_work(i_my_batch)%size, 1

        i_full_points_2 = i_full_points_2 + 1

        if (partition_tab(i_full_points_2).gt.0.d0) then

           i_point = i_point+1

           ! get current integration point coordinate
           coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

              !-----------shanghui add for moving grids----------------
               current_atom = batches(i_my_batch) % points(i_index) % index_atom
               atom_point(i_point) = current_atom
              !-----------shanghui end add for moving grids----------------


           if(n_periodic > 0)then
              call map_to_center_cell(coord_current(1:3) )
           end if

           ! compute atom-centered coordinates of current integration point,
           ! as viewed from all atoms
           call tab_atom_centered_coords_p0 &
                ( coord_current,  &
                dist_tab_sq(1,i_point),  &
                dir_tab(1,1,i_point), &
                n_centers_integrals, centers_basis_integrals )



           ! determine which basis functions are relevant at current integration point,
           ! and tabulate their indices

           ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
           if (.not.prune_basis_once) then
              call prune_basis_p2 &
                   ( dist_tab_sq(1,i_point), &
                   n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
           endif

        end if
     enddo

     if (prune_basis_once) then
        n_compute_c = batches_work(i_my_batch)%batch_n_compute
        i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
     end if

     ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
     call collect_batch_centers_p2 &
     ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, &
       n_batch_centers, batch_center &
     )

     n_points = i_point


     ! Perform actual integration if more than 0 basis functions
     ! are actually relevant on the present angular shell ...
     if (n_compute_c.gt.0) then

        n_rel_points = 0
        i_point = 0

        ! loop over one batch of integration points
        do i_index = 1, batches_work(i_my_batch)%size, 1

           ! Increment the (global) counter for the grid, to access storage arrays
           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1

              coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)!SAG


              if (flag_rel.eq.REL_zora.or. flag_rel==REL_KOLNING_HARMON) then

                 call tab_global_geometry_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dir_tab(1,1,i_point), &
                      dist_tab_full, &
                      i_r_full, &
                      dir_tab_full_norm, &
                      n_centers_integrals,  centers_basis_integrals)

              end if

              ! for all integrations
              partition(i_point) = partition_tab(i_full_points)
  
              do i_atom=1,n_atoms 
              partition_deriv_local(1:3,i_atom,i_point) =    &  
              partition_deriv(1:3,i_atom,i_full_points)
              enddo

              energy_partition(i_point) = partition_tab(i_full_points)

              ! for vectorized xc
              do i_spin = 1, n_spin, 1
                 local_rho(i_spin,i_point) = rho(i_spin,i_full_points)
              enddo


              if (use_gga) then
                 do i_spin = 1, n_spin, 1
                    do i_coord = 1,3,1
                       local_rho_gradient(i_coord,i_spin,i_point) = &
                            rho_gradient(i_coord,i_spin,i_full_points)
                    enddo
                 enddo
              end if
  
              do i_coord=1,3,1             
              local_v_hartree_gradient(i_coord,1:n_atoms,i_point)=      &
              v_hartree_gradient_std(i_coord,1:n_atoms,i_full_points)

              local_first_order_potential_moving_grid(i_coord,1:n_atoms,i_point)=   &
              first_order_potential_moving_grid_std(i_coord,1:n_atoms,i_full_points)

              dir_tab_backup(i_coord,1:n_atoms,i_point)=   &  
              dir_tab(i_coord,1:n_atoms,i_point)  
              enddo
              dist_tab_backup(1:n_atoms,i_point)=dsqrt(dist_tab_sq(1:n_atoms,i_point))
              dist_tab_sq_backup(1:n_atoms,i_point)=dist_tab_sq(1:n_atoms,i_point)


              n_compute_atoms = 0
              n_compute_fns = 0

              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
              ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
              ! without any copying and without doing any unnecessary operations.
              ! The price is that the interface is no longer explicit in terms of physical
              ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
              !write(use_unit,*) dir_tab(:,1,i_point)
              call prune_radial_basis_p2 &
                   ( n_max_compute_atoms, n_max_compute_fns_ham, &
                     dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
                     n_compute_atoms, atom_index, atom_index_inv, &
                     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                     i_atom_fns, spline_array_start, spline_array_end, &
                     n_centers_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                     fn_atom, n_zero_compute, zero_index_point &
                    )
              !write(use_unit,*) dir_tab(:,1,i_point)
              !write(use_unit,*) '-shanghui test dir_tab-------------'

              ! Tabulate distances, unit vectors, and inverse logarithmic grid units
              ! for all atoms which are actually relevant
              call tab_local_geometry_p2 &
                   ( n_compute_atoms, atom_index, &
                     dist_tab(1,i_point), i_r )

              ! compute trigonometric functions of spherical coordinate angles
              ! of current integration point, viewed from all atoms
              call tab_trigonom_p0 &
                   ( n_compute_atoms, dir_tab(1,1,i_point), trigonom_tab )


!shanghui-gga: we alwayse need gradients
!              if ((use_gga) .or. (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON) ) then
                 ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                 call tab_gradient_ylm_p0  &
                      ( trigonom_tab(1,1), basis_l_max,   &
                      l_ylm_max, n_compute_atoms, atom_index,  &
                      ylm_tab(1,1),   &
                      dylm_dtheta_tab(1,1),   &
                      scaled_dylm_dphi_tab(1,1)  )

!              else
!                ! tabulate distance and Ylm's w.r.t. other atoms
!                call tab_wave_ylm_p0 &
!                   ( n_compute_atoms, atom_index,  &
!                   trigonom_tab, basis_l_max,  &
!                   l_ylm_max, ylm_tab )
!              end if

              ! Now evaluate radial functions
              ! from the previously stored compressed spline arrays
              call evaluate_radial_functions_p0  &
                   ( spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_wave_ordered, radial_wave,  &
                   .false. , n_compute_c, n_max_compute_fns_ham )

              ! tabulate total wave function value for each basis function
              call evaluate_waves_p2  &
                   ( n_compute_c, n_compute_atoms, n_compute_fns, &
                     l_ylm_max, ylm_tab, one_over_dist_tab,   &
                     radial_wave, wave(1,i_point), &
                     rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point &
                   )

              ! in the remaining part of the subroutine, some decisions (scalar
              !  relativity) depend on the potential; must therefore evaluate the
              ! potential and derived quantities right here

              ! Local exchange-correlation parts of the potential are evaluated
              ! right here, to avoid having to store them separately elsewhere.
              ! For large systems, savings are significant
              if ((use_hartree_fock).and.(first_integration)) then
                 call evaluate_xc &
                      ( rho(1,i_full_points), &
                      rho_gradient(1,1,i_full_points), &
                      en_density_xc(i_point), local_xc_derivs(1,i_point), &
                      xc_gradient_deriv(1,1,i_point), local_dVxc_drho(1,i_point), .true., coord_current &
                               )  !SAG
              else
                !---------We need local_dVxc_drho, but evaluate_xc do not have--
                !---------so we need evaluate_xc_shanghui below----------------
                !---------just like integrate_first_order_H-------------------- 
                ! call evaluate_xc  &
                !      ( rho(1,i_full_points),   &
                !      rho_gradient(1,1,i_full_points),  &
                !      en_density_xc(i_point), &
                !      en_density_x, en_density_c, &
                !      local_xc_derivs(1,i_point),  &
                !      xc_gradient_deriv(1,1,i_point), local_dVxc_drho(1,i_point), &
                !      coord_current &
                !      )           !SAG

                 call evaluate_xc_shanghui  &
                      ( rho(1,i_full_points),   &
                      rho_gradient(1,1,i_full_points),  &
                      en_density_xc(i_point), &
                      en_density_x, en_density_c, &
                      local_xc_derivs(1,i_point),  &
                      xc_gradient_deriv(1,1,i_point), local_dVxc_drho(1,i_point), &
                      coord_current &
                      )  

              end if


              do i_spin = 1, n_spin, 1
                 local_potential_parts(i_spin) =   &
                     !--------V_Ze------------------------------
                     !-1.0d0/dsqrt(dist_tab_sq(1,i_point)) &
                    ! -1.0d0/dsqrt(dist_tab_sq(2,i_point))!+& 
                    !--------end VZe---------------------------
                    !---------V_ee-----------------------------
                      hartree_potential(i_full_points)   +   &
                    !  1.0d0/dsqrt(dist_tab_sq(1,i_point)) &
                    ! +1.0d0/dsqrt(dist_tab_sq(2,i_point))!+& 
                    !---------end V_ee-------------------------
                    !---------V_xc-----------------------------
                      local_xc_derivs(i_spin,i_point)
                    !---------end V_xc-------------------------


                 if (use_gga) then
                    sum_of_local_gradients(1:3,i_spin) =   &
                         xc_gradient_deriv(1:3,i_spin,i_point)*4.d0
                 else
                    sum_of_local_gradients(1:3,i_spin) = 0.d0
                 end if


              enddo

              ! Check whether relativistic corrections are needed at the present point.
              ! The check is based entirely on the local parts of the potential - i.e.
              ! in a GGA, the terms due to d(rho*exc)/d(|grad(rho|^2) is not evaluated.
              ! Hopefully this approximation to the full ZORA energy is small.
              if (flag_rel.eq.REL_zora.or. (flag_rel==REL_KOLNING_HARMON)) then

                 ! if we need ZORA, must get the _full_ local geometry in order to
                 ! create the superposition of atomic potentials which is used to estimate
                 ! the potential gradient for ZORA

                 call evaluate_pot_superpos_p0  &
                      (   &
                      i_r_full,   &
                      zora_potential_parts(1),  &
                      n_centers_integrals, centers_basis_integrals )

                 do i_spin = 1, n_spin, 1

                    ! factor 2.d0 required because a factor 1/2 is already included in kinetic_wave later ...
                    zora_operator(i_spin) =  &
                         2.d0 * light_speed_sq /  &
                         ( 2 * light_speed_sq -  &
                         zora_potential_parts(i_spin) )

                 enddo

              end if

!              if ((use_gga) .or. (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON)) then
!----------shanghui-gga: always calculate gradient---------------------
                ! we require the gradient of each basis function
                 ! tabulate radial derivatives of those radial functions
                 ! which are actually non-zero at current point, using vectorized splines
                 call evaluate_radial_functions_p0  &
                      ( spline_array_start, spline_array_end,  &
                      n_compute_atoms, n_compute_fns,   &
                      dist_tab(1,i_point), i_r,  &
                      atom_index, i_basis_fns_inv,  &
                      basis_deriv_ordered,   &
                      radial_wave_deriv(1), .true.,  &
                      n_compute_c, n_max_compute_fns_ham )

               if(.false.) then  !shanghui add for hessian_pulay, we always need hessian

                 ! and finally, assemble the actual gradients
                 call evaluate_wave_gradient_p2  &
                    ( n_compute_c, n_compute_atoms, n_compute_fns, &
                      one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                      l_ylm_max, ylm_tab,  &
                      dylm_dtheta_tab,  &
                      scaled_dylm_dphi_tab,  &
                      radial_wave,  &
                      radial_wave_deriv,  &
                      gradient_basis_wave_npoints(1:n_compute_c,1:3,i_point),  &
                      rad_index, wave_index, l_index, l_count, fn_atom, &
                      n_zero_compute, zero_index_point &
                    )
 
               else
 
                  call evaluate_radial_functions_deriv_p0 &
                      (spline_array_start, spline_array_end, &
                       n_compute_atoms, n_compute_fns,  &
                       dist_tab(1,i_point), i_r(1), &
                       atom_index, i_basis_fns_inv, &
                       basis_deriv_ordered,  &
                       radial_wave_2nd_deriv(1),  &
                       .true., n_max_compute_fns_dens &
                      )
        
                  call evaluate_cartesians &
                       (dir_tab(1,1,i_point), basis_l_max, l_wave_max,  &
                        atom_index, n_compute_atoms, cartesians(1,0,1) )
        
                 call evaluate_cartesian_hessian_and_gradient_terms & 
                       (basis_l_max, l_wave_max, cartesians(1,0,1), &
                        atom_index, n_compute_atoms,  &
                        ylm_tab(1,1),  &
                        sum_gradient(1,1,1),  &
                        sum_hessian(1,1,1))


                 call evaluate_wave_gradient_cartesian_p1 &
                     (dist_tab(1,i_point), i_r(1),  &
                     dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                     n_compute_c, i_basis, atom_index_inv,  &
                     i_basis_fns_inv, radial_wave(1),  &
                     radial_wave_deriv(1), &
                     ylm_tab(1,1),  &
                     sum_gradient(1,1,1),  &
                     gradient_basis_wave_npoints(1:n_compute_c,1:3,i_point) & 
                     ,n_compute_atoms)

                 call evaluate_wave_hessian_cartesian_p0 &
                       (dist_tab(1,i_point), i_r(1),  &
                       dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                       n_compute_c, i_basis, atom_index_inv,  &
                       i_basis_fns_inv, radial_wave(1),  &
                       radial_wave_deriv(1),  &
                       radial_wave_2nd_deriv(1),  &
                       ylm_tab(1,1),  &
                       sum_gradient(1,1,1),  &
                       sum_hessian(1,1,1), &
                       hessian_basis_wave(1,1,i_point), n_compute_atoms)

               endif ! calc_hessian_pulay



!---------------shanghui calclate kinetic-gradient_phi---------------------
              ! First, obtain radial kinetic energy terms from vectorized splines
                  call evaluate_radial_functions_p0  &
                       ( spline_array_start, spline_array_end,  &
                       n_compute_atoms, n_compute_fns,   &
                       dist_tab(1,i_point), i_r,  &
                       atom_index, i_basis_fns_inv,  &
                       basis_kinetic_ordered, kinetic_wave(1),  &
                       .false., n_compute_c, n_max_compute_fns_ham )
     
                      !kinetic_wave =  kinetic_wave *0.5

                  call evaluate_radial_functions_deriv_p0 &
                       (spline_array_start, spline_array_end, &
                       n_compute_atoms, n_compute_fns,  &
                       dist_tab(1,i_point), i_r(1), &
                       atom_index, i_basis_fns_inv, &
                       basis_kinetic_ordered,  &
                       kinetic_wave_deriv(1),  &
                       .true., n_max_compute_fns_dens &
                       )

                      !kinetic_wave_deriv =  kinetic_wave_deriv *0.5

                  call evaluate_wave_gradient_cartesian_p1 &
                      (dist_tab(1,i_point), i_r(1), &
                      dir_tab(1,1,i_point), index_lm,  l_wave_max, &
                      n_compute_c, i_basis, atom_index_inv, &
                      i_basis_fns_inv, kinetic_wave(1), &
                      kinetic_wave_deriv(1), &
                      ylm_tab(1,1), &
                      sum_gradient(1,1,1), &
                      kinetic_gradient_basis_wave(1:n_compute_c,1:3,i_point), & 
                      n_compute_atoms)
                 
!---------------shanghui end calclate kinetic_gradient_phi---------------------


                do i_spin=1,n_spin
                   do i_compute=1,n_compute_c
                   H_times_gradient_psi(i_compute,1:3,i_point,i_spin)= & 
                   kinetic_gradient_basis_wave(i_compute,1:3,i_point)+& 
                   local_potential_parts(i_spin)*  & 
                   gradient_basis_wave_npoints(i_compute,1:3,i_point)
                   enddo
                enddo 


              ! Now, evaluate vector of components H*phi(i,r)
              ! Local potential parts first; in the case of GGA,
              ! the real gradient parts are added further below
              !               if ( (flag_rel/=1)) then
              ! Non-relativistic treatment - simply evaluate
              ! H*phi(i,r) all in one

!!!              ! First, obtain radial kinetic energy terms from vectorized splines
!!!              call evaluate_radial_functions_p0  &
!!!                   ( spline_array_start, spline_array_end,  &
!!!                   n_compute_atoms, n_compute_fns,   &
!!!                   dist_tab(1,i_point), i_r,  &
!!!                   atom_index, i_basis_fns_inv,  &
!!!                   basis_kinetic_ordered, kinetic_wave(1),  &
!!!                   .false., n_compute_c, n_max_compute_fns_ham )


              do i_spin = 1, n_spin, 1
                 call evaluate_H_psi_p2  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab,  &
                   radial_wave, H_times_psi(1, i_point, i_spin),  &
                   local_potential_parts(i_spin),  &
                   kinetic_wave, zora_operator(i_spin), &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point &
                 )


              enddo

              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0


              if ((flag_rel.eq.REL_zora).or. flag_rel==REL_KOLNING_HARMON) then

                 ! Scalar relativistic treatment.
                 ! count number of "truly" relativistic points for ZORA treatment
                 ! of kinetic energy ...

                 do i_spin = 1, n_spin, 1

                    zora_operator(i_spin) =  &
                         light_speed_sq /  &
                         (2 * light_speed_sq -  &
                         zora_potential_parts(i_spin))**2

                    call  add_zora_gradient_part_p0(   &
                         sum_of_local_gradients(1,i_spin),  &
                         i_r_full,  &
                         dir_tab_full_norm,   &
                         dist_tab_full,  &
                         zora_operator(i_spin), &
                         n_centers_integrals, centers_basis_integrals )

                 end do

                 do i_spin = 1, n_spin, 1

                    ! Evaluate difference of scalar relativistic kinetic energy operator for the
                    ! true potential and the superposition of free atom potentials separately, and
                    ! only for all relativistic points in shell. Here, use partially
                    ! integrated version, leading to a vector:
                    ! zora_operator(r)*grad(phi(r,i))

                    zora_operator(i_spin) =  &
                         light_speed_sq *  &
                         (local_potential_parts(i_spin) -  &
                         zora_potential_parts(i_spin))/  &
                         ( 2 * light_speed_sq -  &
                         local_potential_parts(i_spin))/  &
                         ( 2 * light_speed_sq -  &
                         zora_potential_parts(i_spin))

                    call evaluate_zora_vector_p1  &
                         ( zora_operator(i_spin),  &
                         partition_tab(i_full_points),  &
                         gradient_basis_wave(1,1),  &
                         n_compute_c,  &
                         zora_vector1(1, 1, n_rel_points+1, i_spin),  &
                         zora_vector2(1, 1, n_rel_points+1, i_spin), &
                         n_max_compute_ham, t_zora(i_spin)  )

                 enddo

                 if (n_spin.eq.1) then
                   if(t_zora(1)) then
                      n_rel_points = n_rel_points + 1
                   end if
                 else if (n_spin.eq.2) then
                   if(t_zora(1).or.t_zora(2)) then
                      n_rel_points = n_rel_points + 1
                   end if
                 end if

              end if  ! end ZORA preparations

              ! If using a GGA, add the true gradient terms to the Hamiltonian vector
              !if (use_gga .or. (n_rel_points.gt.0)) then
              if (use_gga .or.(flag_rel.eq.REL_zora).or. flag_rel==REL_KOLNING_HARMON) then

                 do i_spin = 1, n_spin, 1
                    call add_gradient_part_to_H_p0  &
                         ( n_compute_c,   &
                         gradient_basis_wave(1,1),  &
                         sum_of_local_gradients(1,i_spin),  &
                         H_times_psi(1, i_point, i_spin) )
                 enddo
              end if
   

           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch


!-------------------------shanghui add for first_H--------------------------- 
             call  evaluate_first_order_rho( &
                   n_points,n_occ_states,       &
                   KS_ev_compute,n_compute_c, i_basis, & 
                   wave,gradient_basis_wave_npoints, &
                   first_order_density_matrix,first_order_DM_rho,local_first_order_rho)

             call  evaluate_first_order_rho_moving_grid( &
                   n_points,n_occ_states,       &
                   KS_ev_compute,n_compute_c, i_basis, & 
                   wave,gradient_basis_wave_npoints, &
                   first_order_density_matrix,local_first_order_rho_moving_grid)
!-------------------------shanghui end add-----------------------------------

!-------------------------shanghui use local_first_order_rho to get gradient_diple_electron----------------
       i_point = 0
       ! loop over one batch of integration points
       do i_index = 1, batches_work(i_my_batch)%size, 1
          ! Increment the (global) counter for the grid, to access storage arrays
          i_full_points_first_order_rho = i_full_points_first_order_rho + 1
          coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

       if (partition_tab(i_full_points_first_order_rho).gt.0.d0) then

           i_point = i_point+1
           do i_atom=1,n_atoms
           do i_coord=1,3

             do j_coord=1,3 

             if(use_partition_deriv) then

              if(j_coord.eq.i_coord) then 
             gradient_dipole(3*i_atom+i_coord-3,j_coord)=gradient_dipole(3*i_atom+i_coord-3,j_coord) &
            -partition_tab(i_full_points_first_order_rho)  &
            *local_first_order_rho_moving_grid(i_coord,i_atom,i_point)*coord_current(j_coord)  & 
            -partition_tab(i_full_points_first_order_rho)*local_rho(1,i_point)    &
            -partition_deriv(i_coord,i_atom,i_full_points_first_order_rho)  &
            *local_rho(1,i_point)*coord_current(j_coord) 
              else 
             gradient_dipole(3*i_atom+i_coord-3,j_coord)=gradient_dipole(3*i_atom+i_coord-3,j_coord) &
            -partition_tab(i_full_points_first_order_rho)  &
            *local_first_order_rho_moving_grid(i_coord,i_atom,i_point)*coord_current(j_coord) & 
            -partition_deriv(i_coord,i_atom,i_full_points_first_order_rho)  & 
            *local_rho(1,i_point)*coord_current(j_coord) 
              endif

             else

             gradient_dipole(3*i_atom+i_coord-3,j_coord)=gradient_dipole(3*i_atom+i_coord-3,j_coord) &
            -partition_tab(i_full_points_first_order_rho)  &
            *local_first_order_rho(i_coord,i_atom,i_point)*coord_current(j_coord) 
 
             endif
           
             enddo

           enddo
           enddo
       endif
       enddo
!------------------------shanghui end use local_first_order_rho to get gradient_diple_electron----------------


        !write(use_unit,*)  'use_partition_deriv---------------->',use_partition_deriv

        if(use_partition_deriv) then   

        do i_spin = 1, n_spin, 1

          call  evaluate_hellman_feynman_hessian(hellman_feynman_hessian,        &
                n_points,partition, partition_deriv_local, dist_tab_backup, dir_tab_backup,       &
                local_rho(i_spin,1:n_points),atom_point,local_first_order_rho, &  
                local_first_order_rho_moving_grid)


           ! H_pulay(1)_{uv}:  
            call evaluate_first_order_H_pulay(first_order_H_pulay, &
                 n_points,partition, n_compute_c, i_basis,  &
                 H_times_psi(1,1,i_spin), gradient_basis_wave_npoints)

           !H_pulay(2)_{uv} :
            call evaluate_second_order_H_pulay(second_order_H_pulay,  &
                 partition,dist_tab_backup, dir_tab_backup, &
                 n_points,n_compute_c, i_basis, index_hessian, &
                 wave, gradient_basis_wave_npoints,H_times_psi(1,1,i_spin), & 
                 H_times_gradient_psi(1,1,1,i_spin), hessian_basis_wave,  &  
                 local_first_order_rho,local_v_hartree_gradient, & 
                 local_dVxc_drho(i_spin,1:n_points) , &
                 atom_point, partition_deriv_local, & 
                 local_first_order_rho_moving_grid,local_first_order_potential_moving_grid)

           !S(1)_{uv}
            call evaluate_first_order_S(first_order_S, partition, &
                 n_points, n_compute_c,i_basis,   &
                 wave, gradient_basis_wave_npoints)

           !S(2)_{uv} 
            call evaluate_second_order_S(second_order_S, partition, &
                 n_points,n_compute_c, i_basis, index_hessian, &
                 wave, gradient_basis_wave_npoints,hessian_basis_wave, & 
                 atom_point, partition_deriv_local)

        enddo

        else  
  
        do i_spin = 1, n_spin, 1

          call  evaluate_hellman_feynman_hessian_fixed(hellman_feynman_hessian,        &
                n_points,partition, dist_tab_backup, dir_tab_backup,       &
                local_rho(i_spin,1:n_points),local_first_order_rho)


           ! H_pulay(1)_{uv}:  
            call evaluate_first_order_H_pulay(first_order_H_pulay, &
                 n_points,partition, n_compute_c, i_basis,  &
                 H_times_psi(1,1,i_spin), gradient_basis_wave_npoints)

           !H_pulay(2)_{uv} :
            call evaluate_second_order_H_pulay_fixed(second_order_H_pulay,  &
                 partition,dist_tab_backup, dir_tab_backup, &
                 n_points,n_compute_c, i_basis, index_hessian, &
                 wave, gradient_basis_wave_npoints,H_times_psi(1,1,i_spin), & 
                 H_times_gradient_psi(1,1,1,i_spin), hessian_basis_wave,  &  
                 local_first_order_rho,local_v_hartree_gradient, & 
                 local_dVxc_drho(i_spin,1:n_points))

           !S(1)_{uv}
            call evaluate_first_order_S(first_order_S, partition, &
                 n_points, n_compute_c,i_basis,   &
                 wave, gradient_basis_wave_npoints)

           !S(2)_{uv} 
            call evaluate_second_order_S_fixed(second_order_S, partition, &
                 n_points,n_compute_c, i_basis, index_hessian, &
                 wave, gradient_basis_wave_npoints,hessian_basis_wave) 

        enddo

        endif ! use_partition_deriv

     else

       i_full_points = i_full_points + batches_work(i_my_batch)%size
       i_full_points_first_order_rho = i_full_points_first_order_rho + batches_work(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

     if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  end do ! end loop over batches


   
    !-------shanghui begin parallel------
    if(.not. use_local_index) call sync_integrate_hellman_feynman_hessian( hellman_feynman_hessian )

    if(.not. use_local_index) call sync_integrate_first_order_S( first_order_S )
    if(.not. use_local_index) call sync_integrate_first_order_H( first_order_H_pulay )
    if(.not. use_local_index) call sync_integrate_second_order_S( second_order_S )
    if(.not. use_local_index) call sync_integrate_second_order_H_pulay(second_order_H_pulay )
    !-------shanghui end parallel------
    




!!!debug       write(use_unit,*) '----------------------H_pulay---------------------------------------------------------'
!!!debug       write(use_unit,*) '************shanghui begain first_order_H_pulay(x 1)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_H_pulay(1,1,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_H_pulay****************'
!!!debug       write(use_unit,*) '************shanghui begain first_order_H_pulay(x 2)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_H_pulay(1,2,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_H_pulay****************'
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '************shanghui begain second_order_H_pulay(x 1 x 1)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (second_order_H_pulay(1,1,1,1,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end second_order_H_pulay****************'
!!!debug       write(use_unit,*) '************shanghui begain second_order_H_pulay(x 2 x 2)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (second_order_H_pulay(1,2,1,2,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end second_order_H_pulay****************'
!!!debug
!!!debug
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '----------------------S---------------------------------------------------------'
!!!debug       write(use_unit,*) '************shanghui begain first_order_S(x 1)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_S(1,1,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_S****************'
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '************shanghui begain first_order_S(x 2)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_S(1,2,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_S****************'
!!!debug       write(use_unit,*) '************shanghui begain second_order_S(x 1 x 1)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (second_order_S(1,1,1,1,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end second_order_S****************'
!!!debug
!!!debug       write(use_unit,*) '************shanghui begain second_order_S(x 2 x 2)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (second_order_S(1,2,1,2,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end second_order_S****************'
!!!debug
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '----------------------DM---------------------------------------------------------'
!!!debug       write(use_unit,*) '************shanghui begain DM****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (density_matrix(i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end DM****************'
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '************shanghui begain first_order_DM(1 1)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_density_matrix(1,1,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_DM****************'
!!!debug       write(use_unit,*) '************shanghui begain first_order_DM(1 2)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_density_matrix(1,2,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_DM****************'
!!!debug
!!!debug
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '----------------------EDM---------------------------------------------------------'
!!!debug       write(use_unit,*) '************shanghui begain EDM****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (energy_density_matrix(i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end EDM****************'
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) '************shanghui begain first_order_EDM(x 1)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_energy_density_matrix(1,1,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_EDM****************'
!!!debug       write(use_unit,*) '************shanghui begain first_order_EDM(x 2)****************'
!!!debug       do i_compute=1,n_basis
!!!debug        write(use_unit,'(6f15.9)') (first_order_energy_density_matrix(1,2,i_compute,j_basis),j_basis=1,n_basis )
!!!debug       enddo
!!!debug       write(use_unit,*) '************shanghui end first_order_EDM****************'
!!!debug
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) 'shanghui test second_order_H_pulay:'
!!!debug       write(use_unit,*) ' x 1 x 1: u=1,v=2', second_order_H_pulay(1,1,1,1,1,2)
!!!debug       write(use_unit,*) ' x 1 x 2: u=1,v=2', second_order_H_pulay(1,1,1,2,1,2)
!!!debug       write(use_unit,*) ' x 2 x 1: u=1,v=2', second_order_H_pulay(1,2,1,1,1,2)
!!!debug       write(use_unit,*) ' x 2 x 2: u=1,v=2', second_order_H_pulay(1,2,1,2,1,2)
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) ' x 1 x 1: u=1,v=1', second_order_H_pulay(1,1,1,1,1,1)
!!!debug       write(use_unit,*) ' x 1 x 2: u=1,v=1', second_order_H_pulay(1,1,1,2,1,1)
!!!debug       write(use_unit,*) ' x 2 x 1: u=1,v=1', second_order_H_pulay(1,2,1,1,1,1)
!!!debug       write(use_unit,*) ' x 2 x 2: u=1,v=1', second_order_H_pulay(1,2,1,2,1,1)
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) ' x 1 x 1: u=2,v=2', second_order_H_pulay(1,1,1,1,2,2)
!!!debug       write(use_unit,*) ' x 1 x 2: u=2,v=2', second_order_H_pulay(1,1,1,2,2,2)
!!!debug       write(use_unit,*) ' x 2 x 1: u=2,v=2', second_order_H_pulay(1,2,1,1,2,2)
!!!debug       write(use_unit,*) ' x 2 x 2: u=2,v=2', second_order_H_pulay(1,2,1,2,2,2)
!!!debug       write(use_unit,*) ''
!!!debug------------------end test for second_order_S-----------------------------
!!!debug       
!!!debug
!!!debug         write(use_unit,*) ''
!!!debug         write(use_unit,*) 'shanghui test hellman_feynman_hessian:TOTAL'
!!!debug         write(use_unit,*) ' x 1 x 1',hellman_feynman_hessian(1,1,1,1)+ionic_hessian(1,1,1,1)
!!!debug         write(use_unit,*) ' x 1 x 2',hellman_feynman_hessian(1,1,1,2)+ionic_hessian(1,1,1,2)
!!!debug         write(use_unit,*) ' x 2 x 1',hellman_feynman_hessian(1,2,1,1)+ionic_hessian(1,2,1,1)
!!!debug         write(use_unit,*) ' x 2 x 2',hellman_feynman_hessian(1,2,1,2)+ionic_hessian(1,2,1,2)
!!!debug         write(use_unit,*) ' y 1 y 1', hellman_feynman_hessian(2,1,2,1)+ &
!!!debug                ionic_hessian(2,1,2,1)
!!!debug
!!!debug
!!!debug         write(use_unit,*) ''
!!!debug         write(use_unit,*) 'shanghui test hellman_feynman_hessian==<1>===ionic'
!!!debug         write(use_unit,*) ' x 1 x 1',ionic_hessian(1,1,1,1)
!!!debug         write(use_unit,*) ' x 1 x 2',ionic_hessian(1,1,1,2)
!!!debug         write(use_unit,*) ' x 2 x 1',ionic_hessian(1,2,1,1)
!!!debug         write(use_unit,*) ' x 2 x 2',ionic_hessian(1,2,1,2)
!!!debug         write(use_unit,*) ' y 1 y 1',ionic_hessian(2,1,2,1)
!!!debug
!!!debug         write(use_unit,*) ''
!!!debug         write(use_unit,*) 'shanghui test hellman_feynman_hessian:===<2>===VZe'
!!!debug         write(use_unit,*) ' x 1 x 1',hellman_feynman_hessian(1,1,1,1)
!!!debug         write(use_unit,*) ' x 1 x 2',hellman_feynman_hessian(1,1,1,2)
!!!debug         write(use_unit,*) ' x 2 x 1',hellman_feynman_hessian(1,2,1,1)
!!!debug         write(use_unit,*) ' x 2 x 2',hellman_feynman_hessian(1,2,1,2)
!!!debug         write(use_unit,*) ' y 1 y 1', hellman_feynman_hessian(2,1,2,1)

!------------------test for hellman_feyman_hessian-----------------------------


       !---------shanghui now from all the matrix we have to get Pulay_Hessian: 

           !input denstiy_matrix                   = DM  
           !input first_order_denstiy_matrix       = DM(1)
           !input energy_density_matrix(u,v)       = EDM 
           !input first_order_energy_density_matrix(u,v,natom*3) = EDM(1)

           !calc first_order_S(u,v,natom*3)          = S(1) 
           !calc second_order_S(u,v,natom*3,natom*3) = S(2)
           !calc first_order_H_pulay(u,v,natom*3)    = H_pulay(1)
           !calc second_order_H_pulay(u,v,natom*3,natom*3)= H_pulay(2)

!----------------------begein Pulay hessian test---------------------------------
         ! -----Pulay_force-----------------
         ! sum_{uv}{ DM_{uv}*H_pulay(1)_{uv} }  
         ! sum_{uv}{ EDM_{uv}*S(1)_{uv} }
          
          do i_atom=1,n_atoms
          do i_coord=1,3 
             pulay_force(i_coord,i_atom)=0.0d0
                do i_basis_1=1,n_basis
                do j_basis=1,n_basis
             pulay_force(i_coord,i_atom)= & 
             pulay_force(i_coord,i_atom)+ &
             density_matrix(i_basis_1,j_basis)*  & 
             first_order_H_pulay(i_coord,i_atom,i_basis_1,j_basis)-   & 
             energy_density_matrix(i_basis_1,j_basis)*  & 
             first_order_S(i_coord,i_atom,i_basis_1,j_basis) 
                enddo
                enddo 
           enddo
           enddo

!!!debug       write(use_unit,*) 'shanghui test pulay_force:'
!!!debug       write(use_unit,*) ' x 1:', pulay_force(1,1)
!!!debug       write(use_unit,*) ' x 2:', pulay_force(1,2)
!!!debug       write(use_unit,*) ''
!----------------------end Pulay hessian test---------------------------------
         
          do i_atom=1,n_atoms
          do i_coord=1,3 
             do j_atom=1,n_atoms
             do j_coord=1,3

             pulay_hessian(i_coord,i_atom,j_coord,j_atom)=0.0d0

                do i_basis_1=1,n_basis
                do j_basis=1,n_basis

         !-------so we have pulay_hessian--->(1) 
         ! sum_{uv}{ dm(1)_{uv}*h_pulay(1)_{uv} }  
         ! -sum_{uv}{ edm(1)_{uv}*s(1)_{uv} }
            pulay_hessian(i_coord,i_atom,j_coord,j_atom)= & 
            pulay_hessian(i_coord,i_atom,j_coord,j_atom)+ &
             first_order_density_matrix(j_coord,j_atom,i_basis_1,j_basis)*  & 
             first_order_H_pulay(i_coord,i_atom,i_basis_1,j_basis)-   & 
             first_order_energy_density_matrix(j_coord,j_atom,i_basis_1,j_basis)*  & 
             first_order_S(i_coord,i_atom,i_basis_1,j_basis)+  & 
         !-------and the pulay_hessian   --->(2) 
         ! sum_{uv}{ dm_{uv}*h_pulay(2)_{uv} }  
         ! -sum_{uv}{ edm_{uv}*s(2)_{uv} }
             density_matrix(i_basis_1,j_basis)* & 
             second_order_H_pulay(i_coord,i_atom,j_coord,j_atom,i_basis_1,j_basis)- &
             energy_density_matrix(i_basis_1,j_basis)* & 
             second_order_S(i_coord,i_atom,j_coord,j_atom,i_basis_1,j_basis)

 
                enddo
                enddo 
            
             enddo !j_coord
             enddo !j_atom
        enddo !i_coord
        enddo !i_atom

 

!------------------test for Pulay_hessian---------------------------------
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) 'shanghui test Pulay_hessian:'
!!!debug       write(use_unit,*) ' x 1 x 1', pulay_hessian(1,1,1,1)
!!!debug       write(use_unit,*) ' x 1 x 2', pulay_hessian(1,1,1,2)
!!!debug       write(use_unit,*) ' x 2 x 1', pulay_hessian(1,2,1,1)
!!!debug       write(use_unit,*) ' x 2 x 2', pulay_hessian(1,2,1,2)
!!!debug       write(use_unit,*) ' y 1 y 1', pulay_hessian(2,1,2,1)
!------------------end test for Pulay_hessian---------------------------------



!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) 'shanghui test total_hessian:'
!!!debug       write(use_unit,*) ' x 1 x 1', hellman_feynman_hessian(1,1,1,1)+ & 
!!!debug                ionic_hessian(1,1,1,1)+ pulay_hessian(1,1,1,1)
!!!debug       write(use_unit,*) ' x 1 x 2', hellman_feynman_hessian(1,1,1,2)+ &
!!!debug                ionic_hessian(1,1,1,2)+pulay_hessian(1,1,1,2)
!!!debug       write(use_unit,*) ' x 2 x 1', hellman_feynman_hessian(1,2,1,1)+ & 
!!!debug                ionic_hessian(1,2,1,1)+pulay_hessian(1,2,1,1)
!!!debug       write(use_unit,*) ' x 2 x 2', hellman_feynman_hessian(1,2,1,2)+ & 
!!!debug                ionic_hessian(1,2,1,2)+pulay_hessian(1,2,1,2)
!!!debug       write(use_unit,*) ' y 1 y 1', hellman_feynman_hessian(2,1,2,1)+ &
!!!debug                ionic_hessian(2,1,2,1)+ pulay_hessian(2,1,2,1)
!!!debug       write(use_unit,*) ' y 1 y 2', hellman_feynman_hessian(2,1,2,2)+ &
!!!debug                ionic_hessian(2,1,2,2)+ pulay_hessian(2,1,2,2)
!!!debug
!!!debug       write(use_unit,*) ''
!!!debug       write(use_unit,*) 'shanghui test total_hessian:(ev/A*A)'
!!!debug       write(use_unit,*) ' x 1 x 1', (hellman_feynman_hessian(1,1,1,1)+ & 
!!!debug                ionic_hessian(1,1,1,1)+ pulay_hessian(1,1,1,1))*  & 
!!!debug                hartree / (bohr*bohr)
!!!debug       write(use_unit,*) ' x 1 x 2', (hellman_feynman_hessian(1,1,1,2)+ &
!!!debug                ionic_hessian(1,1,1,2)+pulay_hessian(1,1,1,2))* & 
!!!debug                hartree / (bohr*bohr)
!!!debug       write(use_unit,*) ' x 2 x 1', (hellman_feynman_hessian(1,2,1,1)+ & 
!!!debug                ionic_hessian(1,2,1,1)+pulay_hessian(1,2,1,1)) * &
!!!debug                hartree / (bohr*bohr)
!!!debug       write(use_unit,*) ' x 2 x 2', (hellman_feynman_hessian(1,2,1,2)+ & 
!!!debug                ionic_hessian(1,2,1,2)+pulay_hessian(1,2,1,2)) * &
!!!debug                hartree / (bohr*bohr)
!!!debug       write(use_unit,*) ' y 1 y 1', (hellman_feynman_hessian(2,1,2,1)+ &
!!!debug                ionic_hessian(2,1,2,1)+ pulay_hessian(2,1,2,1))* &
!!!debug                hartree / (bohr*bohr)
!!!debug       write(use_unit,*) ' y 1 y 2', (hellman_feynman_hessian(2,1,2,2)+ &
!!!debug                ionic_hessian(2,1,2,2)+ pulay_hessian(2,1,2,2))* &
!!!debug                hartree / (bohr*bohr)

!!!debug     do i_atom = 1 , n_atoms 
!!!debug     do i_coord = 1 , 3
!!!debug
!!!debug        write(use_unit,'(6E20.12)')  & 
!!!debug            (((hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom)+ &
!!!debug              ionic_hessian(i_coord,i_atom,j_coord,j_atom)+ & 
!!!debug              pulay_hessian(i_coord,i_atom,j_coord,j_atom))* &
!!!debug                hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
!!!debug 
!!!debug     enddo
!!!debug     enddo 

      if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'Total_hessian in integrate_hessian.f90:'

      do i_atom = 1 , n_atoms 
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  & 
             (((hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom)+ &
               ionic_hessian(i_coord,i_atom,j_coord,j_atom)+ & 
               pulay_hessian(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo
      enddo 
      endif 


      do i_atom = 1 , n_atoms 
      do i_coord = 1 , 3
        do j_atom=1, n_atoms
        do j_coord =1, 3 
        hessian(3*i_atom+i_coord-3,3*j_atom+j_coord-3)= & 
        ( hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom)+ &
          ionic_hessian(i_coord,i_atom,j_coord,j_atom)+ &
          pulay_hessian(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr)     
        enddo 
        enddo 

      enddo
      enddo 


!--------------------shanghui make asr 2013.10.21----------------------
      if(use_ASR.and.myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'Make Acoustic sum rule in integrate_hessian.f90:'

      do i_coord = 1 , 3
      do j_coord =1, 3

        do i_atom = 1 , n_atoms
        hessian(3*i_atom+i_coord-3,3*i_atom+j_coord-3)=0.0d0

        do j_atom=1, n_atoms

        if(j_atom.ne.i_atom) then
        hessian(3*i_atom+i_coord-3,3*i_atom+j_coord-3)=  &
        hessian(3*i_atom+i_coord-3,3*i_atom+j_coord-3)-  &
        hessian(3*i_atom+i_coord-3,3*j_atom+j_coord-3)
        endif

        enddo
        enddo

      enddo
      enddo

     endif

!--------------------shanghui end make asr 2013.10.21----------------------



!--------------------shanghui use rho to calculate  gradient_dipole for IR intensity-----------------
!    gradient_dipole(1:3*n_atoms,1:3) = 0.0d0
!    i_full_points = 0
! 
!      do i_my_batch = 1, n_my_batches_work
! 
!        if(get_batch_weights) time_start = mpi_wtime()
! 
!        ! loop over one batch
!        do i_index = 1, batches_work(i_my_batch)%size, 1
! 
!          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
!          i_full_points = i_full_points + 1
! 
!          ! Only execute if partition_tab is .gt. zero, else
!          ! we can run into 1/0 !!!
!          if (partition_tab(i_full_points).gt.0.d0) then
! 
!            ! get current integration point coordinate
!            coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)
! 
!           do i_atom =1,n_atoms
!           do i_coord = 1, 3
! 
!           gradient_dipole(3*i_atom+i_coord-3,1:3)=gradient_dipole(3*i_atom+i_coord-3,1:3) & 
!           +partition_tab(i_full_points)  & 
!           *first_order_rho(i_coord,i_atom,i_full_points)*coord_current(1:3)
! 
!           enddo 
!           enddo
! 
!          end if ! end if (partition_tab.gt.0.d0)
!        end do  ! end loop over points in a batch
!        if(get_batch_weights) batch_times(i_my_batch) = batch_times(i_my_batch) + mpi_wtime() - time_start
!      end do ! end loop over batches

  !-------shanghui begin parallel------
  if(.not. use_local_index) call sync_integrate_gradient_dipole( gradient_dipole )
  !-------shanghui end parallel--------

   !-------------add the gradient_dipole_nuclue ---------------------------------
       do i_atom =1,n_atoms 
       do i_coord =1 ,3 
          do j_coord=1,3

            if(j_coord.eq.i_coord) then
             gradient_dipole(3*i_atom+i_coord-3,j_coord)=gradient_dipole(3*i_atom+i_coord-3,j_coord) &
            +species_z(species(i_atom))               
            endif

          enddo 
       enddo
       enddo 
   !------------end add the gradient_dipole_nuclue---------------------------------



  if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'grdient_dipole in integrate_hessian.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3
         write(use_unit,'(3E20.12)')  gradient_dipole(3*i_atom+i_coord-3,1:3)
      enddo
      enddo  
  endif



!--------------------shanghui end use rho calculate gradient_dipole for IR intensity-----------------



  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  if(time_all>time_work*1.3 .and. .not.use_load_balancing) &
    info_str = trim(info_str) // ' => Consider using load balancing!'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  !     synchronise the hamiltonian
!  if(.not. use_local_index) call sync_integrate_hamiltonian( hamiltonian )
  

  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)



  if(allocated( zora_vector2         )) deallocate( zora_vector2         )
  if(allocated( zora_vector1         )) deallocate( zora_vector1         )
  if(allocated( i_r_full             )) deallocate( i_r_full             )
  if(allocated( dir_tab_full_norm    )) deallocate( dir_tab_full_norm    )
  if(allocated( dist_tab_full        )) deallocate( dist_tab_full        )
  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( wave                 )) deallocate( wave                 )
  if(allocated( local_first_order_rho      )) deallocate( local_first_order_rho      )
  if(allocated( local_first_order_rho_moving_grid      )) deallocate( local_first_order_rho_moving_grid      )
  if(allocated( first_order_DM_rho      )) deallocate( first_order_DM_rho      )
  if(allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if(allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( H_times_psi          )) deallocate( H_times_psi          )
  if(allocated( hamiltonian_shell    )) deallocate( hamiltonian_shell    )
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if(allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if(allocated( gradient_basis_wave  )) deallocate( gradient_basis_wave  )
  if(allocated( gradient_basis_wave_npoints)) deallocate( gradient_basis_wave_npoints)
  if(allocated( local_rho_gradient   )) deallocate( local_rho_gradient   )
  if(allocated( local_rho            )) deallocate( local_rho            )
  if(allocated( xc_gradient_deriv    )) deallocate( xc_gradient_deriv    )
  if(allocated( local_xc_derivs      )) deallocate( local_xc_derivs      )
  if(allocated( local_dVxc_drho      )) deallocate( local_dVxc_drho      )
  if(allocated( local_v_hartree_gradient  )) deallocate( local_v_hartree_gradient     )
  if(allocated( local_first_order_potential_moving_grid )) deallocate( local_first_order_potential_moving_grid    )
  if(allocated( en_density_xc        )) deallocate( en_density_xc        )
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )
  if(allocated( dir_tab_backup       )) deallocate( dir_tab_backup       )
  if(allocated( dist_tab_backup      )) deallocate( dist_tab_backup      )
  if(allocated( dist_tab_sq_backup      )) deallocate( dist_tab_sq_backup      )

  if(use_batch_permutation > 0) then
    deallocate(rho)
    deallocate(hartree_potential)
    deallocate(rho_gradient) ! always allocated
    deallocate(ins_idx)
  endif

  if(get_batch_weights) deallocate(batch_times)



end subroutine integrate_hessian
!******
