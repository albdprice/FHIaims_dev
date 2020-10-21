  Subroutine CC_3d_init_intl()

  Use timing
  Use CC_3d
  use localorb_io, only: use_unit
  use geometry, only: species

  Implicit None

  Integer :: j_tmp,errnum
  Double precision :: CC_time_start,CC_time_end,CC_clock_start,CC_clock_end

  ! Allocate all vectors
  ! For RI coefficients
  Allocate(CC_RI_L(CC_mem_bas(CC_mpi_did+1),CC_n_state, &
                   CC_n_state,n_k_points,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_RI_L in CC')

  CC_RI_L = 0.0D0

  Allocate(CC_RI_R(CC_mem_bas(CC_mpi_did+1),CC_n_state, &
                   CC_n_state,n_k_points,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_RI_R in CC')

  CC_RI_R = 0.0D0

  Call get_timestamps(CC_time_start, CC_clock_start)
  ! Prepare RI-V coefficients
  Call CC_3d_get_RI_coeff()
  Call get_timestamps(CC_time_end, CC_clock_end)

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Three-center(RI-V) coefficients are ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating RI coefficients', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
    write(use_unit,"(2x,A)") 'Allocating integral tensors...'
  end if

  ! Allocate the integral vectors
  ! For (ai|kc) (k,c,kk,i,a,ki,ka)
  Allocate(CC_intl_aikc(CC_n_occ,CC_n_vir,n_k_points, &
                        CC_n_occ,CC_mem_aa_d(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_aikc in CC')

  CC_intl_aikc = 0.0d0
 
  ! For (ki|ac) (k,c,kk,i,a,ki,ka)
  Allocate(CC_intl_kiac(CC_n_occ,CC_n_vir,n_k_points, &
                        CC_n_occ,CC_mem_aa_d(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kiac in CC')

  CC_intl_kiac = 0.0d0
 
  ! For (ac|kd) (k,d,kk,c,a,kc,ka)
  Allocate(CC_intl_ackd(CC_n_occ,CC_n_vir,n_k_points, &
                        CC_n_vir,CC_mem_aa_d(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_ackd in CC')

  CC_intl_ackd = 0.0d0
 
  ! For (ac|kd)_nd ---- (k,c<d,a,kc,ka,ka+kk,1)
  !     (kc|ad)_nd ---- (k,c<d,a,kc,kk,ka+kk,2)
  Allocate(CC_intl_ackd_nd(CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3),CC_n_vir, &
                           n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1),2),stat=errnum)
  Call check_allocation(errnum,'CC_intl_ackd_nd in CC')

  CC_intl_ackd_nd = 0.0D0

  ! For (ac|kd)_d ---- (c=d,kc,a,k,ka=kk)
  Allocate(CC_intl_ackd_d(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_vir, &
                          n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_ackd_d in CC')

  CC_intl_ackd_d = 0.0D0

  ! For (li|kc) (k,c,kk,i,l,ki,kl)
  Allocate(CC_intl_likc(CC_n_occ,CC_n_vir,n_k_points, &
                        CC_n_occ,CC_mem_ii_D(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_likc in CC')

  CC_intl_likc = 0.0d0

  ! For (ki|lc) (k,c,kk,i,l,ki,kl)
  Allocate(CC_intl_kilc(CC_n_occ,CC_n_vir,n_k_points, &
                        CC_n_occ,CC_mem_ii_D(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kilc in CC')

  CC_intl_kilc = 0.0d0

  ! For (ki|lc) (k<=l,kk,i,c,ki,kk+kl,1)
  ! For (ki|lc) (k<=l,kk,j,c,ki,kk+kl,2)
  Allocate(CC_intl_kilc_T(CC_mem_ij_D(CC_mpi_did+1,2),n_k_points,CC_n_occ,CC_n_vir, &
                          n_k_points,CC_mem_k1(CC_mpi_gid+1),2),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kilc_T in CC')

  CC_intl_kilc_T = 0.0d0

  ! For (ki|lj) (k<=l,kk,i,j,ki,kk+kl)
  Allocate(CC_intl_kilj(CC_mem_ij_D(CC_mpi_did+1,2),n_k_points,CC_n_occ,CC_n_occ, &
                        n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kilj_I in CC')

  CC_intl_kilj = 0.0d0

  ! For (ai|bc) (c,b,kb,i,a,ki,ka)
  Allocate(CC_intl_aibc(CC_n_vir,CC_n_vir,n_k_points, &
                        CC_n_occ,CC_mem_aa_d(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_aibc in CC')

  CC_intl_aibc = 0.0d0

  ! For (ai|kj) (j,k,kk,i,a,ki,ka)
  Allocate(CC_intl_aikj(CC_n_occ,CC_n_occ,n_k_points, &
                        CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                        CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_aikj in CC')

  CC_intl_aikj = 0.0d0

  if (CC_sv_strategy.eq.1) then

    ! For (kc|ld)_A (l,d,kl,k,c,kk,kc)
    Allocate(CC_intl_kcld_A(CC_n_occ,CC_n_vir,n_k_points, &
                            CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                            n_k_points),stat=errnum)
    Call check_allocation(errnum,'CC_intl_kcld_A in CC')

    CC_intl_kcld_A = 0.0d0

  else if (CC_sv_strategy.eq.2) then
 
    ! For (kc|ld) (l,d,kl,k,c,kk,kc)
    Allocate(CC_intl_kcld(CC_n_occ,CC_n_vir,n_k_points, &
                          CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                          CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
    Call check_allocation(errnum,'CC_intl_kcld in CC')

    CC_intl_kcld = 0.0d0

  end if

  ! For (ia|jb)_nd (i,j,a<b,ki,ka,ka+kb)
  Allocate(CC_intl_iajb_nd(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                           n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_iajb_nd in CC')
  
  CC_intl_iajb_nd = 0.0d0

  ! For (ia|jb)_d (i,j,a=b,ki,ka,ka+kb)
  Allocate(CC_intl_iajb_d(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                          n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_iajb_d in CC')

  CC_intl_iajb_d = 0.0D0

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensors allocated successfully.'
  end if

  ! Calculate integral tensors
  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_kiac()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ki|ac) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ki|ac)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_aikc()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ai|kc) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ai|kc)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_ackd()

  Call CC_3d_intl_ackd_nd()

  Call CC_3d_intl_ackd_d()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ac|kd) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ac|kd)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_kilc()

  Call CC_3d_intl_kilc_T()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ki|lc) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ki|lc)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_likc()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (li|kc) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (li|kc)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_kilj()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ki|lj) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ki|lj)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_aibc()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ai|bc) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ai|bc)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_aikj()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ai|kj) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ai|kj)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  Call CC_3d_intl_iajb()

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (ia|jb) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (ia|jb)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  Call get_timestamps(CC_time_start, CC_clock_start)

  if (CC_sv_strategy.eq.1) then

    Call CC_3d_intl_kcld_A()

  else if (CC_sv_strategy.eq.2) then

    Call CC_3d_intl_kcld()

  end if

  Call get_timestamps(CC_time_end, CC_clock_end)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral tensor (kc|ld) is ready.'
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                 'Time of calculating integral (kc|ld)', &
                 CC_time_end - CC_time_start,'s (cpu)', &
                 CC_clock_end - CC_clock_start, 's (wall clock)'
  end if

  if (CC_abcd_sv_flag) then
    ! (ac|kd) (c<d,kc,a,b,ka,ka+kb) 
    Allocate(CC_intl_acbd_nd(CC_mem_ab_D(CC_mpi_did+1,3),n_k_points,CC_n_vir, &
                             CC_n_vir,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
    Call check_allocation(errnum,'CC_intl_acbd_nd in CC')

    CC_intl_acbd_nd = 0.0d0

    ! (ac|kd) (c=d,kc,a,b,ka,ka+kb)
    Allocate(CC_intl_acbd_d(CC_mem_aa_D(CC_mpi_did+1),n_k_points,CC_n_vir, &
                            CC_n_vir,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
    Call check_allocation(errnum,'CC_intl_acbd_d in CC')

    CC_intl_acbd_d = 0.0d0

    Call get_timestamps(CC_time_start, CC_clock_start)

    Call CC_3d_intl_acbd()

    Call get_timestamps(CC_time_end, CC_clock_end)
    if (myid.eq.0) then
      write(use_unit,"(2x,A)") 'Integral tensor (ac|bd) is ready.'
      write(use_unit,"(2x,A,F14.3,A,F14.3,A)") &
                   'Time of calculating integral (ac|bd)', &
                   CC_time_end - CC_time_start,'s (cpu)', &
                   CC_clock_end - CC_clock_start, 's (wall clock)'
    end if

  end if

  End Subroutine CC_3d_init_intl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_get_RI_coeff()

  Use basis
  Use prodbas
  Use hartree_fock
  Use CC_3d
  Use mpi_tasks
  use geometry, only: species
  use pbc_lists, only: Cbasis_to_atom

  Implicit None

  Integer :: i_k_point,id_root,i_k_point_local,k1,k2,k3,k4,k_12,k_34
  Integer :: k_start,k_end,k_run,n_k1
  Integer :: i_k_pattern,i_k_start,i_run,i_go,i_RI,n_RI,i_k1,i_k2
  Integer :: n_omp,oid,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  Integer :: errnum,g_recv,g_index,i_1,i_2
  Double complex , dimension(:,:,:,:) , allocatable :: lvl_recip1_tmp, &
                                                       coulomb_matr_recip_tmp, &
                                                       KS_egnvec_tmp, &
                                                       lvl_comm_tmp, coulomb_comm_tmp, &
                                                       KS_egnvec_comm_tmp, &
                                                       RI_tmp_R, RI_tmp_L

  Double complex , dimension(:,:) , allocatable :: lvl_tricoeff_sum_all,lvl_nn_tmp, &
                                                   lvl_nn_tmp_2,lvl_nn_tmp_3,mul_tmp


  Integer (kind = 8) :: s_tmp
  Integer :: n_row,n_col
  Integer :: i_row,i_col,j_tmp
  Integer :: i_basis_1,i_basis_2,i_atom_1,i_atom_2,i_species_1,i_species_2
  Integer :: bboff_1,bboff_2,n_spbb_1,n_spbb_2

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  Allocate(lvl_recip1_tmp(max_n_basbas_sp,n_basis,n_basis,n_k_points),stat=errnum)
  Call check_allocation(errnum,'lvl_recip1_tmp in CC')

  lvl_recip1_tmp = 0.0D0

  Allocate(lvl_comm_tmp(max_n_basbas_sp,n_basis,n_basis,1),stat=errnum)
  Call check_allocation(errnum,'lvl_comm_tmp in CC')

  lvl_comm_tmp = 0.0D0

  Allocate(coulomb_matr_recip_tmp(CC_mem_bas(CC_mpi_did+1),n_basbas, &
                                  n_k1,1),stat=errnum)
  Call check_allocation(errnum,'coulomb_matr_recip_tmp in CC')

  coulomb_matr_recip_tmp = 0.0D0

  Allocate(KS_egnvec_tmp(n_basis,n_states,n_k_points,1),stat=errnum)
  Call check_allocation(errnum,'KS_egnvec_tmp in CC')

  KS_egnvec_tmp = 0.0D0

  Allocate(KS_egnvec_comm_tmp(n_basis,n_states,1,1),stat=errnum)
  Call check_allocation(errnum,'KS_egnvec_comm_tmp in CC')
  
  KS_egnvec_comm_tmp = 0.0D0

  if (flag_KS_eigenfunc_conjg) then
    if (allocated(KS_eigenvector_complex)) then
      KS_eigenvector_complex = conjg(KS_eigenvector_complex)
    end if
  end if

  do i_k_point = 1, n_k_points

    id_root = Mod(i_k_point,n_tasks)
    i_k_point_local = Int((i_k_point - 1) / n_tasks) + 1

    if (myid.eq.id_root) then

      lvl_comm_tmp(:,:,:,1) = lvl_tricoeff_recip1(:,:,:,i_k_point_local)

      if (real_eigenvectors) then
        KS_egnvec_comm_tmp(:,:,1,1) = KS_eigenvector(:,:,1,i_k_point_local)
      else
        KS_egnvec_comm_tmp(:,:,1,1) = KS_eigenvector_complex(:,:,1,i_k_point_local)
      end if
    end if

    s_tmp = Int(size(lvl_comm_tmp),8)
    Call CC_mpi_complex_bcast(s_tmp, lvl_comm_tmp, id_root, MPI_COMM_WORLD)

    s_tmp = Int(size(KS_egnvec_comm_tmp),8)
    Call CC_mpi_complex_bcast(s_tmp, KS_egnvec_comm_tmp, id_root, MPI_COMM_WORLD)
        
    lvl_recip1_tmp(:,:,:,i_k_point) = lvl_comm_tmp(:,:,:,1)

    KS_egnvec_tmp(:,:,i_k_point,1) = KS_egnvec_comm_tmp(:,:,1,1)

  end do

  Deallocate(lvl_comm_tmp)
  Deallocate(KS_egnvec_comm_tmp)
 
  g_recv = 0
  g_index = 1

  do i_k_point = 1, n_k_points

    id_root = Mod(i_k_point,n_tasks)
    i_k_point_local = Int((i_k_point - 1) / n_tasks) + 1

    j_tmp = CC_index_k1(g_index) - 1 + CC_mem_k1(g_index)
    if (i_k_point.gt.j_tmp) then
      g_index = g_index + 1
      g_recv = g_index - 1
    end if

    k1 = i_k_point - CC_index_k1(CC_mpi_gid+1) + 1

    !write(80+myid,*) 'i_k_point',i_k_point
    !write(80+myid,*) 'g_index',g_index
    !write(80+myid,*) 'g_recv',g_recv
    !write(80+myid,*) 'k1',k1
    !write(80+myid,*) 'id_root',id_root

    if (myid.eq.id_root) then
      j_tmp = g_recv
      i_run = 1
      do while (j_tmp.le.n_tasks-1)
 
        i_1 = CC_index_bas(i_run)
        i_2 = i_1 - 1 + CC_mem_bas(i_run)

        s_tmp = Int(CC_mem_bas(i_run) * n_basbas,8)
        if (id_root.ne.j_tmp) then
          Call CC_mpi_complex_send(s_tmp,coulomb_matr_recip(i_1:i_2,:,i_k_point_local), &
                                   j_tmp,i_k_point,MPI_COMM_WORLD)
        else
          coulomb_matr_recip_tmp(:,:,k1,1) = coulomb_matr_recip(i_1:i_2,:,i_k_point_local)
        end if

        i_run = i_run + 1
        j_tmp = j_tmp + CC_n_domain
      end do
    else if (CC_mpi_gid.eq.g_recv) then

      s_tmp = Int(CC_mem_bas(CC_mpi_did+1) * n_basbas,8)
      Call CC_mpi_complex_recv(s_tmp,coulomb_matr_recip_tmp(:,:,k1,1), &
                               id_root,i_k_point,MPI_COMM_WORLD)
    end if

  end do

  Allocate(RI_tmp_R(CC_mem_bas(CC_mpi_did+1),n_basis,n_basis,1),stat=errnum)
  Call check_allocation(errnum,'RI_tmp_R in CC')

  Allocate(RI_tmp_L(CC_mem_bas(CC_mpi_did+1),n_basis,n_basis,1),stat=errnum)
  Call check_allocation(errnum,'RI_tmp_L in CC')

  Allocate(mul_tmp(CC_mem_bas(CC_mpi_did+1),1),stat=errnum)
  Call check_allocation(errnum,'mul_tmp in CC')

  Allocate(lvl_nn_tmp(n_basis,n_basis),stat=errnum)
  Call check_allocation(errnum,'lvl_nn_tmp in CC')

  Allocate(lvl_nn_tmp_2(n_basis,n_states),stat=errnum)
  Call check_allocation(errnum,'lvl_nn_tmp in CC')

  Allocate(lvl_nn_tmp_3(n_states,n_states),stat=errnum)
  Call check_allocation(errnum,'lvl_nn_tmp in CC')

  Allocate(lvl_tricoeff_sum_all(n_basbas,1),stat=errnum)
  Call check_allocation(errnum,'lvl_tricoeff_sum_all in CC')

  do k1 = 1, n_k_points

    do k_run = 1, n_k1

      ! k_12 = k2 - k1
      k_12 = k_start - 1 + k_run

      Call CC_3d_k_plus(k2,k1,k_12)

      do i_basis_1 = 1, n_basis
        do i_basis_2 = 1, n_basis
  
          i_atom_1 = Cbasis_to_atom(i_basis_1)
          i_species_1 = species(i_atom_1)
          bboff_1 = atom2basbas_off(i_atom_1)
          n_spbb_1 = sp2n_basbas_sp(i_species_1)
  
          i_atom_2 = Cbasis_to_atom(i_basis_2)
          i_species_2 = species(i_atom_2)
          bboff_2 = atom2basbas_off(i_atom_2)
          n_spbb_2 = sp2n_basbas_sp(i_species_2)
  
          ! Calculate CC_RI_L
          lvl_tricoeff_sum_all(:,1) = 0.0D0
  
          lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,1) = &
            lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,1) + &
            conjg(lvl_recip1_tmp(1:n_spbb_1, i_basis_1, i_basis_2,k1)) + &
            lvl_tricoeff_recip2(1:n_spbb_1, i_basis_2, i_basis_1)
  
          lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,1) = &
            lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,1) + &
            lvl_recip1_tmp(1:n_spbb_2, i_basis_2, i_basis_1,k2) + &
            lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)
 
 
          i_RI = CC_index_bas(CC_mpi_did+1) - 1
  
          do i_go = 1, CC_mem_bas(CC_mpi_did+1)
            i_RI = i_RI + 1
            RI_tmp_L(i_go,i_basis_2,i_basis_1,1) = &
                     lvl_tricoeff_sum_all(i_RI,1)
          end do
 
          ! Calculate CC_RI_R
          Call zgemm('N', 'N', CC_mem_bas(CC_mpi_did+1), 1, n_basbas, (1.d0,0.d0), &
                     coulomb_matr_recip_tmp(:,:,k_run,1), CC_mem_bas(CC_mpi_did+1), &
                     lvl_tricoeff_sum_all, n_basbas, (0.d0,0.d0), &
                     mul_tmp, CC_mem_bas(CC_mpi_did+1))
 
          RI_tmp_R(:,i_basis_2,i_basis_1,1) = mul_tmp(:,1)

        end do
      end do

!      print*,'RI_tmp_L'
!      print*,RI_tmp_L
!      print*,'RI_tmp_R'
!      print*,RI_tmp_R

      do i_go = 1, CC_mem_bas(CC_mpi_did+1)
  
        lvl_nn_tmp = RI_tmp_L(i_go,:,:,1)
        Call zgemm('N', 'N', n_basis, n_states, n_basis, (1.d0,0.d0), &
                   lvl_nn_tmp, n_basis, KS_egnvec_tmp(:,:,k2,1), n_basis, &
                   (0.d0,0.d0), lvl_nn_tmp_2, n_basis)
  
        Call zgemm('C', 'N', n_states, n_states, n_basis, (1.d0,0.d0), &
                   KS_egnvec_tmp(:,:,k1,1), n_basis, lvl_nn_tmp_2, n_basis, &
                   (0.d0,0.d0), lvl_nn_tmp_3, n_states)
 
!        print*,'lvl_nn_tmp_3'
!        print*,lvl_nn_tmp_3
!        print*,'conjg(lvl_nn_tmp_3'
!        print*,conjg(lvl_nn_tmp_3)

        CC_RI_L(i_go,:,:,k1,k2) = lvl_nn_tmp_3(CC_n_fc+1:n_states,CC_n_fc+1:n_states)
  
        lvl_nn_tmp = RI_tmp_R(i_go,:,:,1)
        Call zgemm('N', 'N', n_basis, n_states, n_basis, (1.d0,0.d0), &
                   lvl_nn_tmp, n_basis, KS_egnvec_tmp(:,:,k2,1), n_basis, &
                   (0.d0,0.d0), lvl_nn_tmp_2, n_basis)
  
        Call zgemm('C', 'N', n_states, n_states, n_basis, (1.d0,0.d0), &
                   KS_egnvec_tmp(:,:,k1,1), n_basis, lvl_nn_tmp_2, n_basis, &
                   (0.d0,0.d0), lvl_nn_tmp_3, n_states)
  
        CC_RI_R(i_go,:,:,k1,k2) = lvl_nn_tmp_3(CC_n_fc+1:n_states,CC_n_fc+1:n_states) &
                                        / dble(n_k_points)
  
      end do

    end do
  end do

  Deallocate(lvl_tricoeff_sum_all)
  Deallocate(lvl_nn_tmp,lvl_nn_tmp_2,lvl_nn_tmp_3)
  Deallocate(lvl_recip1_tmp,coulomb_matr_recip_tmp,KS_egnvec_tmp)

  j_tmp = CC_mem_bas(CC_mpi_did+1) * CC_n_state**2
  s_tmp = Int(j_tmp,8) * Int(n_k_points**2,8)

  Call CC_mpi_complex_allreduce(s_tmp, CC_RI_L, CC_mpi_comm_group)
  Call CC_mpi_complex_allreduce(s_tmp, CC_RI_R, CC_mpi_comm_group)

!  print*,'RI_L'
!  print*,CC_RI_L
!  print*,'RI_R'
!  print*,CC_RI_R

  End Subroutine CC_3d_get_RI_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_kiac()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp,ccc

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (ki|ac) (k{k3}i{k2}|a{k1}c{k4})
  ! CC_intl_kiac (k,c,k3,i,a,k2,k1)
  CC_intl_kiac = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_occ
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_vir 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_occ,n_dom,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do  k_run = 1, n_k1

      k1 = k1 + 1
      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k3,k2,k1,k4,4)

          intl_tmp1(:,:,:) = CC_RI_L(:,1:CC_n_occ,1:CC_n_occ,k3,k2)
          intl_tmp2(:,:,:) = CC_RI_R(:,a_start:a_end,CC_n_occ+1:CC_n_occ+CC_n_vir,k1,k4)

          shp1(1) = nnn
          shp1(2) = i_tmp
         
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_occ
          shp2i(3) = n_dom
          shp2i(4) = CC_n_vir

          intl_tmp3 = reshape(intl_C,shp2i)

          ! (k{k3}i{k2}|a{k1}c{k4}) --> (k,c,k3,i,a,k2,k1)
          do ccc = 1, CC_n_vir
            intl_rlt(:,ccc,k3,:,:,k2,k_run) = intl_tmp3(:,:,:,ccc)
          end do

        end do
      end do

    end do

    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_kiac = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_kiac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_aikc()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (ai|kc) (a{k1}i{k2}|k{k3}c{k4})
  ! CC_intl_aikc(k,c,k3,i,a,k2,k1)

  CC_intl_aikc = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_vir 
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_occ 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_vir,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do k_run = 1, n_k1

      k1 = k1 + 1

      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1(:,:,:) = CC_RI_R(:,1:CC_n_occ,CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)
          intl_tmp2(:,:,:) = CC_RI_L(:,a_start:a_end,1:CC_n_occ,k1,k2)

          shp1(1) = nnn
          shp1(2) = i_tmp
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_vir
          shp2i(3) = n_dom
          shp2i(4) = CC_n_occ
          intl_tmp3 = reshape(intl_C,shp2i)

          do iii = 1, CC_n_occ
            do aaa = 1, n_dom
              intl_rlt(:,:,k3,iii,aaa,k2,k_run) = intl_tmp3(:,:,aaa,iii)
            end do
          end do

        end do
      end do

    end do
   
    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_aikc = intl_rlt 
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_aikc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_ackd()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp,ccc

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (ac|kd) (a{k1}c{k2}|k{k3}d{k4})
  ! CC_intl_ackd (k,d,k3,c,a,k2,k1)
  CC_intl_ackd = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_vir,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_vir 
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_vir 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_vir,n_dom,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do k_run = 1, n_k1

      k1 = k1 + 1

      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1(:,:,:) = CC_RI_R(:,1:CC_n_occ,CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)
          intl_tmp2(:,:,:) = CC_RI_L(:,a_start:a_end,CC_n_occ+1:CC_n_occ+CC_n_vir,k1,k2)

          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_vir
          shp2i(3) = n_dom
          shp2i(4) = CC_n_vir

          intl_tmp3 = reshape(intl_C,shp2i)

          ! (a{k1}c{k2}|k{k3}d{k4}) --> (k,d,k3,c,a,k2,k1)
          do ccc = 1, CC_n_vir
            intl_rlt(:,:,k3,ccc,:,k2,k_run) = intl_tmp3(:,:,:,ccc)
          end do

        end do
      end do

    end do

    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_ackd = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_ackd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_ackd_nd()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,kkk,aaa,bbb,ccc,ddd,code_ia,code_jb,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end,code_cd,c_run
  Integer :: target_id,i_use,i_task

  ! For integrals (ac|kd) (a{k1}c{k2}|k{k3}d{k4}) & (k{k1}c{k2}|a{k3}d{k4})
  ! For (ac|kd)_nd ---- (k,c<d,a,kc,ka,ka+kk,1)
  !     (kc|ad)_nd ---- (k,c<d,a,kc,kk,ka+kk,2)
  CC_intl_ackd_nd = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_ab_D(i_task,3)
    a_start = CC_index_ab_D(i_task,3)
    a_end = a_start - 1 + n_dom

    Allocate(intl_rlt(CC_n_occ,n_dom,CC_n_vir,n_k_points, &
                      n_k_points,n_k1,2),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_vir,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do  k_run = 1, n_k1

      k_13 = k_13 + 1
      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          code_cd = a_start - 1
          do c_run = 1, n_dom

            code_cd = code_cd + 1

            Call CC_3d_decode(code_cd,ccc,ddd,CC_n_vir,3)

            ccc = ccc + CC_n_occ
            ddd = ddd + CC_n_occ

            intl_A = CC_RI_L(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ccc,k1,k2)
            intl_B = CC_RI_R(:,1:CC_n_occ,ddd,k3,k4)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_vir,CC_n_occ,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_vir)

            ! For (ac|kd)_nd ---- (k,c<d,a,kc,ka,ka+kk,1)
            do kkk = 1, CC_n_occ
              intl_rlt(kkk,c_run,:,k2,k1,k_run,1) = intl_C(:,kkk)
            end do

            intl_A = CC_RI_L(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ddd,k3,k4)
            intl_B = CC_RI_R(:,1:CC_n_occ,ccc,k1,k2)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_vir,CC_n_occ,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_vir)

            !     (kc|ad)_nd ---- (k,c<d,a,kc,kk,ka+kk,2)
            do kkk = 1, CC_n_occ
              intl_rlt(kkk,c_run,:,k2,k1,k_run,2) = intl_C(:,kkk)
            end do

          end do

        end do
      end do
    end do

    s_tmp = Int(n_dom,8) * Int(CC_n_vir,8) * Int(CC_n_occ,8) &
          * Int(n_k1*2,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_ackd_nd = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  End Subroutine CC_3d_intl_ackd_nd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_ackd_d()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,kkk,lll,ccc,ddd,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end,code_cd,c_run
  Integer :: target_id,i_use,i_task

  ! For integrals (ac|kd) (a{k1}c{k2}|k{k3}d{k4}) & (k{k1}c{k2}|a{k3}d{k4})
  ! For (ac|kd)_d ---- (k,c=d,a,kc,ka,ka+kk)
  CC_intl_ackd_d = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task)
    a_end = a_start - 1 + n_dom

    Allocate(intl_rlt(CC_n_occ,n_dom,CC_n_vir,n_k_points, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_vir,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do  k_run = 1, n_k1

      k_13 = k_13 + 1
      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          code_cd = a_start - 1
          do c_run = 1, n_dom

            code_cd = code_cd + 1

            ccc = code_cd + CC_n_occ
            ddd = ccc

            intl_A = CC_RI_L(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ccc,k1,k2)
            intl_B = CC_RI_R(:,1:CC_n_occ,ddd,k3,k4)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_vir,CC_n_occ,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_vir)

            ! For (ac|kd)_d ---- (k,c=d,a,kc,ka,ka+kk,1)
            do kkk = 1, CC_n_occ
              intl_rlt(kkk,c_run,:,k2,k1,k_run) = intl_C(:,kkk)
            end do

          end do

        end do
      end do
    end do

    s_tmp = Int(n_dom,8) * Int(CC_n_vir,8) * Int(CC_n_occ,8) &
          * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_ackd_d = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  End Subroutine CC_3d_intl_ackd_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_kilc()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,ccc,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,i_start,i_end
  Integer :: target_id,i_use,i_task

  ! For integrals (ki|lc) (k{k3}i{k2}|l{k1}c{k4})
  ! CC_intl_kilc(k,c,k3,i,l,k2,k1)
  CC_intl_kilc = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_ii_D(i_task)
    i_start = CC_index_ii_D(i_task)
    i_end = i_start - 1 + CC_mem_ii_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_occ
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_vir 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_occ,n_dom,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do  k_run = 1, n_k1

      k1 = k1 + 1
      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k3,k2,k1,k4,4)

          intl_tmp1(:,:,:) = CC_RI_L(:,1:CC_n_occ,1:CC_n_occ,k3,k2)
          intl_tmp2(:,:,:) = CC_RI_R(:,i_start:i_end,CC_n_occ+1:CC_n_occ+CC_n_vir,k1,k4)

          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_occ
          shp2i(3) = n_dom
          shp2i(4) = CC_n_vir

          intl_tmp3 = reshape(intl_C,shp2i)

          ! (k{k3}i{k2}|l{k1}c{k4}) --> CC_intl_kilc(k,c,k3,i,l,k1,k2)
          do ccc = 1, CC_n_vir
            intl_rlt(:,ccc,k3,:,:,k2,k_run) = intl_tmp3(:,:,:,ccc)
          end do

        end do
      end do

    end do

    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_kilc = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_kilc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_kilc_T()



  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,kl_start,kl_end,code_kl,kkk,lll,kl_run
  Integer :: target_id,i_use,i_task

  ! For integrals (ki|lc) (k{k1}i{k2}|l{k3}c{k4})
  ! CC_intl_kilc(ki|lc) (k<=l,kk,i,c,ki,kk+kl,1)
  ! CC_intl_kilc(ki|lc) (k<=l,kk,j,c,ki,kk+kl,2)
  CC_intl_kilc_T = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_ij_D(i_task,2)
    kl_start = CC_index_ij_D(i_task,2)
    kl_end = kl_start - 1 + n_dom

    Allocate(intl_rlt(n_dom,n_k_points,CC_n_occ,CC_n_vir, &
                      n_k_points,n_k1,2),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do  k_run = 1, n_k1

      k_13 = k_13 + 1

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          code_kl = kl_start - 1
          do kl_run = 1, n_dom

            code_kl = code_kl + 1
            Call CC_3d_decode(code_kl,kkk,lll,CC_n_occ,2)

            intl_A = CC_RI_L(:,kkk,1:CC_n_occ,k1,k2)
            intl_B = CC_RI_R(:,lll,CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_occ,CC_n_vir,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_occ)

            intl_rlt(kl_run,k1,:,:,k2,k_run,1) = intl_C

            intl_A = CC_RI_L(:,lll,1:CC_n_occ,k3,k4)
            intl_B = CC_RI_R(:,kkk,CC_n_occ+1:CC_n_occ+CC_n_vir,k1,k2)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_occ,CC_n_vir,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_occ)

            intl_rlt(kl_run,k1,:,:,k2,k_run,2) = intl_C

          end do

        end do
      end do
    end do

    s_tmp = Int(n_dom,8) * Int(CC_n_occ*CC_n_vir,8) * Int(n_k1*2,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_kilc_T = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  End Subroutine CC_3d_intl_kilc_T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_likc()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,lll,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,i_start,i_end
  Integer :: target_id,i_use,i_task

  ! For integrals (li|kc) (l{k1}i{k2}|k{k3}c{k4})
  ! CC_intl_likc(k,c,k3,i,l,k2,k1)
  CC_intl_likc = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_ii_D(i_task)
    i_start = CC_index_ii_D(i_task)
    i_end = i_start - 1 + CC_mem_ii_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_vir 
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_occ
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_vir,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do k_run = 1, n_k1

      k1 = k1 + 1

      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1 = CC_RI_R(:,1:CC_n_occ,CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)
          intl_tmp2 = CC_RI_L(:,i_start:i_end,1:CC_n_occ,k1,k2)

          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_vir
          shp2i(3) = n_dom
          shp2i(4) = CC_n_occ

          intl_tmp3 = reshape(intl_C,shp2i)

          !(l{k1}i{k2}|k{k3}c{k4}) --> CC_intl_likc(k,c,k3,i,l,k2,k1)
          do lll = 1, n_dom
            intl_rlt(:,:,k3,:,lll,k2,k_run) = intl_tmp3(:,:,lll,:)
          end do

        end do
      end do

    end do
   
    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_likc = intl_rlt 
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_likc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_kilj()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1,k_run
  Integer :: nnn,errnum
  Integer :: i_run,code_kl,kkk,lll

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,i_start,i_end
  Integer :: target_id,i_task

  ! For integrals (ki|lj) (k{k1}i{k2}|l{k3}j{k4})
  ! CC_intl_kilj(k<=l,kk,i,j,ki,kk+kl)
  CC_intl_kilj = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_ij_D(i_task,2)
    i_start = CC_index_ij_D(i_task,2)
    i_end = i_start - 1 + n_dom

    Allocate(intl_rlt(n_dom,n_k_points,CC_n_occ,CC_n_occ, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')
    
    Allocate(intl_B(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do  k_run = 1, n_k1

      k_13 = k_13 + 1

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          code_kl = i_start - 1
          do i_run = 1, n_dom

            code_kl = code_kl + 1
            Call CC_3d_decode(code_kl,kkk,lll,CC_n_occ,2)

            intl_A = CC_RI_L(:,kkk,1:CC_n_occ,k1,k2)
            intl_B = CC_RI_R(:,lll,1:CC_n_occ,k3,k4)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_occ,CC_n_occ,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_occ)

            ! CC_intl_kilj(k,l,k1,i<=j,k2,k1+k3)
            intl_rlt(i_run,k1,:,:,k2,k_run) = intl_C

          end do
        end do
      end do

    end do

    s_tmp = Int(CC_n_occ**2,8) * Int(n_dom,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_kilj = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  End Subroutine CC_3d_intl_kilj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_aibc()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (ai|bc) (a{k1}i{k2}|b{k3}c{k4})
  ! CC_intl_aibc(c,b,k4,i,a,k2,k1)
  CC_intl_aibc = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_vir,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_vir * CC_n_vir 
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_vir,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_occ 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_vir,CC_n_vir,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do  k_run = 1, n_k1

      k1 = k1 + 1
      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1(:,:,:) = CC_RI_R(:,CC_n_occ+1:CC_n_occ+CC_n_vir, &
                                        CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)
          intl_tmp2(:,:,:) = CC_RI_L(:,a_start:a_end,1:CC_n_occ,k1,k2)


          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_vir
          shp2i(2) = CC_n_vir
          shp2i(3) = n_dom
          shp2i(4) = CC_n_occ

          intl_tmp3 = reshape(intl_C,shp2i)

          ! (a{k1}i{k2}|b{k3}c{k4}) --> (c,b,k3,i,a,k1,k2)
          do bbb = 1, CC_n_vir
            do aaa = 1, n_dom
              intl_rlt(:,bbb,k4,:,aaa,k2,k_run) = intl_tmp3(bbb,:,aaa,:)
            end do
          end do

        end do
      end do

    end do

    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_aibc = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_aibc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_aikj()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (ai|kj) (a{k1}i{k2}|k{k3}j{k4})
  ! CC_intl_aikj_W(j,k,k4,i,a,k2,k1)

  CC_intl_aikj = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_occ,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_occ 
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = n_dom * CC_n_occ 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_occ,n_dom,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k1 = k_start - 1
    do k_run = 1, n_k1

      k1 = k1 + 1

      do k3 = 1, n_k_points
        do k2 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1(:,:,:) = CC_RI_R(:,1:CC_n_occ,1:CC_n_occ,k3,k4)
          intl_tmp2(:,:,:) = CC_RI_L(:,a_start:a_end,1:CC_n_occ,k1,k2)

          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_occ
          shp2i(3) = n_dom
          shp2i(4) = CC_n_occ

          intl_tmp3 = reshape(intl_C,shp2i)

          ! (a{k1}i{k2}|k{k3}j{k4}) --> (j,k,k3,i,a,k2,k1)
          do jjj = 1, CC_n_occ
            do aaa = 1, n_dom
              intl_rlt(jjj,:,k4,:,aaa,k2,k_run) = intl_tmp3(:,jjj,aaa,:)
            end do
          end do

        end do
      end do

    end do
   
    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_aikj = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_aikj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_kcld()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (kc|ld) (k{k1}c{k2}|l{k3}d{k4}) 
  ! CC_intl_kcld (l,d,k3,k,c,k1,k2)
  CC_intl_kcld = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_vir 
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = CC_n_occ * n_dom 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,CC_n_occ,n_dom),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_vir,CC_n_occ,n_dom),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k2 = k_start - 1
    do  k_run = 1, n_k1

      k2 = k2 + 1
      do k3 = 1, n_k_points
        do k1 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1(:,:,:) = CC_RI_R(:,1:CC_n_occ,CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)
          intl_tmp2(:,:,:) = CC_RI_L(:,1:CC_n_occ,a_start:a_end,k1,k2)

          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_vir
          shp2i(3) = CC_n_occ
          shp2i(4) = n_dom

          intl_tmp3 = reshape(intl_C,shp2i)

          ! (l,d,k3,k,c,k1,k2)
          intl_rlt(:,:,k3,:,:,k1,k_run) = intl_tmp3(:,:,:,:)

        end do
      end do

    end do

    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_kcld = intl_rlt 
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  End Subroutine CC_3d_intl_kcld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_kcld_A()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ia,code_jb,i_run,i_tmp,j_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  ! For integrals (kc|ld) (k{k1}c{k2}|l{k3}d{k4})
  ! CC_intl_kcld_A (l,d,k3,k,c,k1,k2)
  CC_intl_kcld_A = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + CC_mem_aa_D(i_task)

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    i_tmp = CC_n_occ * CC_n_vir
    Allocate(intl_A(nnn,i_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_tmp1(nnn,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_tmp1 in CC')
    
    j_tmp = CC_n_occ * n_dom 
    Allocate(intl_B(nnn,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_tmp2(nnn,CC_n_occ,n_dom),stat=errnum)
    Call check_allocation(errnum,'intl_tmp2 in CC')

    Allocate(intl_C(i_tmp,j_tmp),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    Allocate(intl_tmp3(CC_n_occ,CC_n_vir,CC_n_occ,n_dom),stat=errnum)
    Call check_allocation(errnum,'intl_tmp3 in CC')

    k2 = k_start - 1
    do  k_run = 1, n_k1

      k2 = k2 + 1
      do k3 = 1, n_k_points
        do k1 = 1, n_k_points

          Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

          intl_tmp1(:,:,:) = CC_RI_R(:,1:CC_n_occ,CC_n_occ+1:CC_n_occ+CC_n_vir,k3,k4)
          intl_tmp2(:,:,:) = CC_RI_L(:,1:CC_n_occ,a_start:a_end,k1,k2)

          shp1(1) = nnn
          shp1(2) = i_tmp
     
          intl_A = reshape(intl_tmp1,shp1)

          shp1(1) = nnn
          shp1(2) = j_tmp
          intl_B = reshape(intl_tmp2,shp1)

          alpha = 1.0D0
          beta = 0.0D0

          Call Zgemm('T','N',i_tmp,j_tmp,nnn,alpha,intl_A,nnn, &
                     intl_B,nnn,beta,intl_C,i_tmp)

          shp2i(1) = CC_n_occ
          shp2i(2) = CC_n_vir
          shp2i(3) = CC_n_occ
          shp2i(4) = n_dom

          intl_tmp3 = reshape(intl_C,shp2i)

          intl_rlt(:,:,k3,:,:,k1,k_run) = intl_tmp3(:,:,:,:)

        end do
      end do

    end do

    s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_kcld_A(:,:,:,:,:,:,k_start:k_end) = intl_rlt
    end if

    Deallocate(intl_rlt,intl_A,intl_B,intl_C,intl_tmp1,intl_tmp2,intl_tmp3)

  end do

  i_tmp = CC_n_occ * CC_n_vir
  j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1)
  s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k_points**3,8)
  Call CC_mpi_complex_allreduce(s_tmp, CC_intl_kcld_A, CC_mpi_comm_group)


  End Subroutine CC_3d_intl_kcld_A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_iajb()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_13,k_24,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: iii,jjj,aaa,bbb,code_ab,i_run,i_tmp,j_tmp,ab_run

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp2i

  Integer (kind = 8) :: s_tmp
  Integer :: n_a,k_start,k_end,a_start,a_end,ab_start,ab_end,n_ab
  Integer :: i_task,target_id

  ! For integrals (ia|jb)_nd (i{k2}a{k1}|j{k4}b{k3})
  ! CC_intl_iajb_nd (i,j,a<b,k2,k1,k1+k3)
  CC_intl_iajb_nd = 0.0D0
  CC_intl_iajb_d = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_ab = CC_mem_ab_D(i_task,3)
    ab_start = CC_index_ab_D(i_task,3)
    ab_end = ab_start - 1 + n_ab

    Allocate(intl_rlt(CC_n_occ,CC_n_occ,n_ab, &
                      n_k_points,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do k_run = 1, n_k1

      k_13 = k_13 + 1

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          code_ab = ab_start - 1
          do ab_run = 1, n_ab

            code_ab = code_ab + 1

            Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)
            aaa = aaa + CC_n_occ
            bbb = bbb + CC_n_occ

            intl_A(:,:) = CC_RI_L(:,1:CC_n_occ,aaa,k2,k1)
            intl_B(:,:) = CC_RI_R(:,1:CC_n_occ,bbb,k4,k3)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_occ,CC_n_occ,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_occ)

            ! (i{k2}a{k1}|j{k4}b{k3}) --> (i,j,a<b,k2,k1,k1+k3)
            intl_rlt(:,:,ab_run,k2,k1,k_run) = intl_C
          end do

        end do
      end do
    end do

    s_tmp = Int(CC_n_occ**2,8) * Int(n_ab,8) * Int(n_k_points**2,8) * Int(n_k1,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_iajb_nd = intl_rlt
    end if

    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

    n_a = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) + CC_n_occ
    a_end = a_start - 1 + n_a

    Allocate(intl_rlt(CC_n_occ,CC_n_occ,n_a, &
                      n_k_points,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do k_run = 1, n_k1

      k_13 = k_13 + 1

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          aaa = a_start - 1
          do ab_run = 1, n_a

            aaa = aaa + 1
            bbb = aaa

            intl_A(:,:) = CC_RI_L(:,1:CC_n_occ,aaa,k2,k1)
            intl_B(:,:) = CC_RI_R(:,1:CC_n_occ,bbb,k4,k3)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_occ,CC_n_occ,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_occ)

            ! (i{k2}a{k1}|j{k4}b{k3}) --> (i,j,a<b,k2,k1,k1+k3)
            intl_rlt(:,:,ab_run,k2,k1,k_run) = intl_C
          end do

        end do
      end do
    end do

    s_tmp = Int(CC_n_occ**2,8) * Int(n_a,8) * Int(n_k_points**2,8) * Int(n_k1,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_iajb_d = intl_rlt
    end if

    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  n_ab = CC_mem_ab_D(CC_mpi_did+1,3)
  ab_start = CC_index_ab_D(CC_mpi_did+1,3)
  ab_end = ab_start - 1 + n_ab

  End Subroutine CC_3d_intl_iajb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_intl_acbd()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_24,k_13,k_pattern,n_k1
  Integer :: nnn,errnum
  Integer :: k_run

  Double complex , dimension(:,:) , allocatable :: intl_A, intl_B, intl_C
  Double complex , dimension(:,:,:) , allocatable :: intl_tmp1, intl_tmp2
  Double complex , dimension(:,:,:,:) , allocatable :: intl_tmp3
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: intl_rlt
  Double complex :: alpha,beta

  Integer :: aaa,bbb,ccc,ddd,i_run,i_tmp,j_tmp,code_cd,c_run

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i

  Integer (kind = 8) :: s_tmp
  Integer :: n_dom,i_dom,k_start,k_end,a_start,a_end
  Integer :: target_id,i_use,i_task

  !print*,'intl_abcd'

  ! For integrals (ac|bd) (a{k1}c{k2}|b{k3}d{k4})
  ! CC_intl_acbd_nd(c<d,k2,a,b,k1,k1+k3)
  CC_intl_acbd_nd = 0.0D0
  CC_intl_acbd_d = 0.0D0
  nnn = CC_mem_bas(CC_mpi_did+1)

  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_ab_D(i_task,3)
    a_start = CC_index_ab_D(i_task,3)
    a_end = a_start - 1 + n_dom

    Allocate(intl_rlt(n_dom,n_k_points,CC_n_vir,CC_n_vir, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_vir,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do  k_run = 1, n_k1

      k_13 = k_13 + 1
      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          code_cd = a_start - 1
          do c_run = 1, n_dom

            code_cd = code_cd + 1

            Call CC_3d_decode(code_cd,ccc,ddd,CC_n_vir,3)

            ccc = ccc + CC_n_occ
            ddd = ddd + CC_n_occ

            intl_A = CC_RI_L(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ccc,k1,k2)
            intl_B = CC_RI_R(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ddd,k3,k4)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_vir,CC_n_vir,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_vir)

            ! CC_intl_acbd(c,d,a,b,k2,k_13)
            intl_rlt(c_run,k2,:,:,k1,k_run) = intl_C

            !Call CC_3d_decode(code_cd,ccc,ddd,CC_n_vir,3)
            !do aaa = 1, CC_n_vir
            !  do bbb = 1, CC_n_vir
            !    print*,aaa,bbb,ccc,ddd,intl_rlt(c_run,k2,aaa,bbb,k1,k_run)
            !  end do
            !end do
          end do

        end do
      end do
    end do

    s_tmp = Int(n_dom,8) * Int(CC_n_vir**2,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_acbd_nd = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  ! For integrals (ac|bd) (a{k1}c{k2}|b{k3}d{k4})
  ! CC_intl_acbd_d(c=d,k2,a,b,k1,k1+k3)
  n_k1 = CC_mem_k1(CC_mpi_gid+1)

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  do i_task = 1, CC_mpi_domain_size

    n_dom = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task)
    a_end = a_start - 1 + n_dom

    Allocate(intl_rlt(n_dom,n_k_points,CC_n_vir,CC_n_vir, &
                      n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')
    intl_rlt = 0.0D0

    Allocate(intl_A(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_A in CC')

    Allocate(intl_B(nnn,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_B in CC')

    Allocate(intl_C(CC_n_vir,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl_C in CC')

    k_13 = k_start - 1
    do  k_run = 1, n_k1

      k_13 = k_13 + 1
      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          code_cd = a_start - 1
          do c_run = 1, n_dom

            code_cd = code_cd + 1

            ccc = code_cd + CC_n_occ
            ddd = ccc

            intl_A = CC_RI_L(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ccc,k1,k2)
            intl_B = CC_RI_R(:,CC_n_occ+1:CC_n_occ+CC_n_vir,ddd,k3,k4)

            alpha = 1.0D0
            beta = 0.0D0

            Call Zgemm('T','N',CC_n_vir,CC_n_vir,nnn,alpha,intl_A,nnn, &
                       intl_B,nnn,beta,intl_C,CC_n_vir)

            ! CC_intl_acbd(c,d,a,b,k2,k_13)
            intl_rlt(c_run,k2,:,:,k1,k_run) = intl_C

            !do aaa = 1, CC_n_vir
            !  do bbb = 1, CC_n_vir
            !    print*,aaa,bbb,code_cd,code_cd,intl_rlt(c_run,k2,aaa,bbb,k1,k_run)
            !  end do
            !end do
 
          end do

        end do
      end do
    end do

    s_tmp = Int(n_dom,8) * Int(CC_n_vir**2,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_intl_acbd_d = intl_rlt
    end if
 
    Deallocate(intl_rlt,intl_A,intl_B,intl_C)

  end do

  End Subroutine CC_3d_intl_acbd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

