
!****s* FHI-aims/diagonalize_auxmat_scalapack_complex
!  NAME
!   power_auxmat_scalapack
!  SYNOPSIS

      subroutine diagonalize_auxmat_scalapack_complex &
                  (n_bb, auxmat, prod_threshold, n_nonsingular, eigenvalues, eigenvectors, sqrtv_eigenvectors)

!  PURPOSE
!   Diagonalize a square matrix within the auxiliary basis.
!   Modified from the subroutine power_auxmat_scalapack_complex.
!   X. Ren, 2017.07.25 
!   This is the scalapack version
!
!   This version works on a 2-D distributed matrix
!
!  USES
       use dimensions
       use runtime_choices
       use prodbas
       use mpi_tasks
       use synchronize_mpi
       use scalapack_wrapper
       use crpa_blacs
       use elpa2_2013
       implicit none

!  ARGUMENTS
      integer, intent(IN) :: n_bb
      complex*16, intent(IN) :: auxmat(lbb_row:ubb_row, lbb_col:ubb_col)
      real*8, intent(IN) :: prod_threshold
      real*8, intent(OUT) :: eigenvalues(n_bb)
      complex*16, intent(OUT) :: eigenvectors(n_bb_row,n_bb_col)
      complex*16, intent(OUT) :: sqrtv_eigenvectors(n_bb_row,n_bb_col)

!
!  INPUTS
!    o auxmat -- real array, e.g. bare coulomb matrix within the auxiliary basis
!    o power -- real number, e.g. -0.5d0
!  OUTPUTS
!    o auxmat -- real array, e.g. the inverse square root of the bare coulomb matrix
!                within the auxiliary basis
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

!      real*8, dimension(:), allocatable ::    eigenvalues
!      complex*16, dimension(:,:), allocatable ::  eigenvectors, auxmat2
      complex*16, dimension(:,:), allocatable :: auxmat2
      real*8  :: ev_sqrt 

!     working array
      integer   n_nonsingular

      integer   info
      integer   mpierr
      real*8 :: diff, diff_all, checkloc, checkglob
      character*150 :: info_str
      complex*16:: zero, one
      integer:: npcol_blacs, mypcol_blacs
      integer:: nprow_blacs, myprow_blacs
!  counter
      integer i_basis_1
      integer lc

!  start to work

      call mpi_comm_size(comm_blacs_row,nprow_blacs,mpierr)
      call mpi_comm_rank(comm_blacs_row,myprow_blacs,mpierr)
      call mpi_comm_size(comm_blacs_col,npcol_blacs,mpierr)
      call mpi_comm_rank(comm_blacs_col,mypcol_blacs,mpierr)

!      if(myid.eq.0 .and. name /= '') then
!       write (use_unit,'(2X,4A)') &
!        "Diagonalizing the ", trim(name), " matrix and ", &
!        "checking for singularities with ScaLapack ."
!       write (use_unit,*)
!      endif


!      allocate(eigenvectors(n_bb_row,n_bb_col))
      allocate(auxmat2(n_bb_row,n_bb_col))
      zero=0.
      one=1.
      eigenvectors=0.

! check if matrix is self-adjoint
      call pztranc( n_bb, n_bb, &
           one, auxmat, 1, 1, bb2desc, &
           zero, eigenvectors, 1, 1, bb2desc )
      checkloc=sum(abs(auxmat-eigenvectors))
      call mpi_allreduce(checkloc, checkglob,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm_blacs,mpierr)
!     XR, 2016.10.14
!     Changing the threshold to a bigger value to avoid a lot of such information printed out
!     in periodic GW calculations. An error of 1.e-9 seems to be harmless. 
!     if((checkglob.gt.5e-11).and.(myid_bl.eq.0)) write (use_unit,'(2X,2A,ES14.4)') &
      if((checkglob.gt.5e-9).and.(myid_bl.eq.0)) write (use_unit,'(2X,2A,ES14.4)') &
           'WARNING: input matrix in diagonalize_auxmat_scalapack_complex does not seem to be self-adjoint!',&
           ' Checksum: ',checkglob

! averaging the upper and lower triangles
      auxmat2(:,:) = 0.5d0* (auxmat(:,:) + eigenvectors(:,:))

! diagonalizing the matrix


      auxmat2(:,:) = - auxmat2(:,:)

      call solve_evp_complex_2stage_2013(n_bb, n_bb, auxmat2, &
           ubound(auxmat2,1), eigenvalues, eigenvectors, &
           ubound(eigenvectors,1), bb_bl_row, comm_blacs_row, &
           comm_blacs_col, comm_blacs)

      eigenvalues(1:n_bb) = -eigenvalues(1:n_bb)

      do i_basis_1 = 1, n_bb, 1
         if(eigenvalues(i_basis_1) < prod_threshold) exit
      enddo

      n_nonsingular = i_basis_1-1
!      n_nonsingular = n_bb

!      if(myid.eq.0 .and. name /= '') then
!         write(use_unit,'(2X,3A,ES14.4,A,ES14.4,A)') &
!         & "Eigenvalues of the ", trim(name), " matrix range from", &
!         & eigenvalues(1), " to", eigenvalues(n_bb), "."
!         write(use_unit,'(2X,A,I8,A,I8,3A)') &
!         & "Using ", n_nonsingular, "   eigenvalues out of rank ", &
!         & n_basbas, "   ", trim(name), " matrix (auxiliary basis)."
!         if (n_nonsingular < n_basbas) then
!            write(use_unit,'(2X,A,ES14.4,A,ES14.4,3A)') &
!            & "Still using eigenvalue ", eigenvalues(n_nonsingular), &
!            & " while cutting ", eigenvalues(n_nonsingular+1), &
!            & " in ", trim(name), " matrix."
!         end if
!         write(use_unit,*)
!      end if

      lc = 0 ! local column counter
      do i_basis_1 = 1, n_nonsingular, 1
         if( MOD((i_basis_1-1)/bb_bl_col,npcol_blacs) .eq. mypcol_blacs) then
            ! column i is on local processor
            lc = lc+1
            ev_sqrt = sqrt(eigenvalues(i_basis_1))
            sqrtv_eigenvectors(:, lc) = &
               eigenvectors(:, lc) * ev_sqrt
         endif
      enddo


!  deallocate arrays

       if(allocated(auxmat2)) then
         deallocate(auxmat2)
       endif

     end subroutine diagonalize_auxmat_scalapack_complex

!******
