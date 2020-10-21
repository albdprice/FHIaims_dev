! slightly simplify exception check from socket
FUNCTION drv_check(info)
    use synchronize_mpi_basic
    IMPLICIT NONE
    LOGICAL drv_check
    INTEGER, INTENT(IN) :: info

    call mp_bcast_int(info)
    drv_check = .false.
    if (info.le.0) drv_check = .true.
END FUNCTION

!MR TODO: Probably good to pass also the converged_cpsf flag, but needs testing.

SUBROUTINE run_driver(converged_scf, enough_walltime)
    use localorb_io
    use mpi_utilities
    use synchronize_mpi_basic
    use analytical_stress
    use geometry
    use runtime_choices
    use relaxation
    use physics
    use species_data 
    use dimensions
    use constants  
    use mpi_tasks
    use potconstrain 
    use heat_flux, only: HF_stress_per_atom
    implicit none
    
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL :: isinit=.true., hasdata=.false., drv_check
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer
    INTEGER socket, nat, siz, info, thisrep
    INTEGER inet, port, ccmd, i, run_counter, replicaid, i_atom
    CHARACTER*1024 :: host  
    REAL *8 :: cellh(3,3), cellih(3,3), vir(3,3), pot
    REAL*8, ALLOCATABLE :: combuf(:)
    ! FlK: Atomic stress
    REAL *8 :: stress_per_atom(3, 3, n_atoms)
    REAL*8 :: sigma(3,3)
    logical :: converged_scf, converged_cpscf, enough_walltime 
    CHARACTER*1024 srvaddress
    CHARACTER*100000 extra_string
    CHARACTER*1024 info_str
    real*8 :: free_energy

    replicaid=0
    thisrep=0
    converged_scf=.false. ! start whole thing without having converged the scf obvs.
    ! this should be an option (or several) in the input
!    srvaddress="localhost:12345"
    srvaddress=trim(pimdwrapper_server)//":"//trim(pimdwrapper_port)
    inet=1    ! 1 if we want a internet socket (remote name or localhost or ip)   0 if we want UNIX socket
    host=srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)//achar(0)    ! host name (or unix socket name)
    read(srvaddress(INDEX(srvaddress,':',back=.true.)+1:),*) port       ! port number (kind of ignored for unix sockets)
    
    IF (myid.eq.0) write(*,*) " @ DRIVER MODE: Connecting to host:port ", &
             trim(srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)), port
    
    IF (srvaddress(1:INDEX(srvaddress,':')-1).eq.('UNIX')) THEN
      inet=0
      host=srvaddress(6:INDEX(srvaddress,':',back=.true.)-1)//achar(0)    
    ENDIF

    IF (myid.eq.0) call open_socket(socket, inet, port, host)          
    ! initialize run counter only because for FHI-aims the initialization of the scf is a bit different
    run_counter=0
    
    driver_loop: DO
      ! In this case the walltime is governed by the wrapper, so FHI-aims doesn't need to know about it - always true.
      enough_walltime=.true.
      ! do communication on master node only...
      if (myid.eq.0) call readbuffer(socket, header, MSGLEN, info)
      if (drv_check(info)) exit
      call mp_bcast_char(header) ! character
      
      if (myid.eq.0) write(*,*) " @ DRIVER MODE: Message from server: ", header
      if (trim(header) == "STATUS") then
         if (myid.eq.0) then  
            if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN, info)
            else if (isinit) then
               call writebuffer(socket,"READY       ",MSGLEN, info)
            else 
               call writebuffer(socket,"NEEDINIT       ",MSGLEN, info)
            endif
         endif
         if (drv_check(info)) exit
      else if (trim(header)=="INIT") then
         if (myid.eq.0) call readbuffer(socket, thisrep, 4, info) ! reading replica info 
         if (drv_check(info)) exit
         call mp_bcast_int(thisrep)
         if (thisrep.ne.replicaid .and. .not. run_counter.eq.0 .and. recalc_dens) then
            need_dens_superpos=.true.
         endif
         if (myid.eq.0) write(*,*) " @ DRIVER MODE: Receiving replica",thisrep ! , replicaid, ipi_need_dens_superpos 
         replicaid = thisrep
         if (myid.eq.0) then
            call readbuffer(socket, nat, 4, info) ! length of parameter string -- ignored at present!
            call readbuffer(socket, parbuffer, nat, info)
        endif
        if (drv_check(info)) exit
        isinit=.true.
      else if (trim(header) == "POSDATA") then              
         if (myid.eq.0) then        
            call readbuffer(socket, cellh, 9*8, info)
            call readbuffer(socket, cellih, 9*8, info)
            call readbuffer(socket, nat, 4, info)
         endif         
         if (drv_check(info)) exit
         call mp_bcast_mat(cellh) ! 2d matrix 
         call mp_bcast_mat(cellih) ! 2d matrix
         call mp_bcast_int(nat) ! integer
         if (.not.allocated(combuf)) allocate(combuf(3*nat))
         if (myid.eq.0) call readbuffer(socket, combuf, nat*3*8, info)
         if (drv_check(info)) exit         
         call mp_bcast_arr(combuf) ! real array
         
         ! converts cell and positions data to internal format
         cellh=transpose(cellh)
         cellih=transpose(cellih)
         coords = RESHAPE(combuf, (/ 3 , nat /) )
         do i=1,n_periodic
           lattice_vector(:,i)=cellh(:,i)
         enddo
         call initialize_bravais_quantities()
         ! This call to map to center cell was necessary for corner cases where FHI-aims decided to make atoms jump to far away cells.
         if (n_periodic.gt.0) then
            do i_atom=1, n_atoms
               call map_to_center_cell(coords(:,i_atom))
            enddo 
         endif
            
         if (myid.eq.0) write(*,*) " @ DRIVER MODE: Received positions "
         
         ! Here is where FHI-aims needs to be called.
          
         ! MR: initialize only if it is the  first run ... 
         if(run_counter.gt.0) then
            call reinitialize_scf ( converged_scf )
         else
            call initialize_scf ( converged_scf )
         endif

         call scf_solver (converged_scf , enough_walltime)

         ! MR: The following block is here to avoid that when the sc_iter_init flag is triggered, i-PI receives buggy forces (all zero)
         if (restarting_scf) then
             if (myid.eq.0) write(*,*) "|   Not updating geometry. Just reinitializing the mixer for convergence" 
             call reinitialize_scf (converged_scf)
             call scf_solver (converged_scf , enough_walltime)
         endif

         ! MR: The IF below is just a failsafe option in case something gets changed in the scf_solver
         ! since the ipi_need_dens_superpos should be set to false after the first iteration in there
         if (need_dens_superpos) need_dens_superpos=.false.
         !
        call clean_force_components(total_forces, forces_lv)
          
         ! MR: If constrain to coordinates, coming from an external potential, is applied, apply here. 
         if (use_potential_constrain) then
            call apply_constrain(free_energy)
            total_energy=total_energy+free_energy
         endif

         ! Starting infrastructure for passing electronic friction
         if(use_friction .and. ipi_ftensor) then
            call calculate_nonadiabatic_friction(converged_scf)
         endif
         ! End infrastructure for passing electronic friction

          ! Start implementing DFPT infrastructure also here
          ! TODO: Probably all these ifs should go to a sort of postprocessing
          ! TODO: Create an ipi flag to pass polarizabilities through the extra infra.
          if(use_DFPT_polarizability) then
             converged_cpscf=.false.
             call cpscf_solver_polarizability( converged_cpscf )
          else if(use_DFPT_dielectric) then
             converged_cpscf=.false.
             call cpscf_solver_dielectric( converged_cpscf )
          endif
         ! End DFPT calcs

         ! converts results to i-pi format
         extra_string = trim(comm_string)  ! additional information as generated by aims 
         comm_string=' ' ! resets comm_string after it has been passed              
         combuf=RESHAPE(total_forces, (/ 3 * nat /) )   ! return force in atomic units

         ! FlK: Save the atomic stress if it was computed
         if (compute_heat_flux) then
            stress_per_atom = HF_stress_per_atom
         else
            stress_per_atom = 0.0d0
         endif

         pot=total_energy + 2d0 * entropy_correction  ! return potential in atomic units         
         if (use_analytical_stress) then
           vir=analytical_stress_tensor ! return virial in atomic units and without the volume scaling
         else if (use_numerical_stress) then
           vir=numerical_stress_tensor
         else
           vir =0.0d0
           if (run_counter.eq.0 .and. myid.eq.0) write(*,*) &
                    " @ DRIVER - WARNING: the stress tensor will not be computed"
         endif         
         vir = -1.0*cell_volume*transpose(vir)  ! scale the virial to be consistent with i-PI

         run_counter=run_counter+1         
         hasdata=.true.
         !! at the end of getfrce
         isinit = .false. ! resets init so that we will get replicaindex again at next step!
         !! Make FHI-aims output the structure it is calculating
         if (myid.eq.0) then
            write (*,*) "------------------------------------------------------------------------"
            write (*,*) "Atomic structure that was used in the preceding time step of the wrapper"
            call output_structure()
            write (*,*) "------------------------------------------------------------------------"
         endif
         !!
      else if (trim(header)=="GETFORCE") then
         if (myid.eq.0) write(*,*) " @ DRIVER MODE: Returning v,forces,stress "
         if (myid.eq.0) then      
            call writebuffer(socket,"FORCEREADY  ",MSGLEN, info)
            call writebuffer(socket,pot,8, info)
            call writebuffer(socket,nat,4, info)            
            call writebuffer(socket,combuf,3*nat*8, info)
            call writebuffer(socket,vir,9*8, info)
          endif
          if (drv_check(info)) exit

          ! MR: it is okay to allow an empty string to be passed to i-PI here since 
          ! i-PI just dumps it and uses it for nothing of note. I may be proven wrong.
          if (myid.eq.0) then
            siz=len(trim(extra_string))
            call writebuffer(socket,siz,4, info)
            call writebuffer(socket,trim(extra_string),siz, info)
         endif
         hasdata=.false.

      ! FlK: Communicate the atomic stress
      else if (trim(header)=="GETSTRESSES") then
         if (myid.eq.0) write(*,*) " @ DRIVER MODE: returning atomic stress "
         if (myid.eq.0) then
            call writebuffer(socket,"STRESSREADY  ",MSGLEN, info)
            call writebuffer(socket,nat,4, info)            
            call writebuffer(socket,stress_per_atom,3*3*nat*8, info)
          endif
          if (drv_check(info)) exit
         hasdata=.false.

      ! FlK: Switch atomic stress computation on and off
      else if (trim(header)=="STRESSES_ON") then
         if (myid.eq.0) then
            write(*,*) " @ DRIVER MODE: switch atomic stress computation on "
         end if
         use_analytical_stress = .true.
         compute_heat_flux     = .true.
      else if (trim(header)=="STRESSES_OFF") then
         if (myid.eq.0) then
            write(*,*) " @ DRIVER MODE: switch atomic stress computation off "
         end if
         use_analytical_stress = .false.
         compute_heat_flux     = .false.
      endif
    ENDDO driver_loop    
    
END SUBROUTINE
  
