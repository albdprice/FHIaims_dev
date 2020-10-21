!****h* FHI-aims/mpe_types
!  NAME
!    mpe_types
!  SYNOPSIS

module mpe_types

!  PURPOSE
!    This module contains all derived types for the MPE continuum solvation
!    module and the corresponding synchronization routines.
!  USES

   use types, only: dp
   use mpi_types, only: type_mpi
   use mpi_headers, only: MPI_ADDRESS_KIND, MPI_INTEGER
   use mpe_constants, only: MPE_CONST

   implicit none

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   private  ! make everything private by default

   ! PUBLIC TYPES
   public   DielectricInterface, InterfacePoint, &
            DielectricContinuum, &
            Basis, SolHarmBasis, &
            BasisCenter, SolHarmBasisCenter, &
            SpVec, SpVecTuple

   ! PUBLIC ROUTINES
   public   initialize_mpe_types, cleanup_mpe_types, &
            SolHarmBasisCenter_size_lm, SolHarmBasisCenter_lbound_lm, SolHarmBasisCenter_ubound_lm, &
            InterfacePoint_vector_mpi_bcast, &
            InterfacePoint_vector_extract, &
            SpVec_mpi_bcast, SpVec_from_dense_vector, &
            SpVecTuple_vector_mpi_bcast, SpVecTuple_vector_from_dense_vector

   ! CONSTANTS
   real(dp), parameter :: ZERO=0.e0_dp, ONE=1.e0_dp

   ! TYPES
   type :: SpVecTuple
         real(dp) :: val
         integer :: ind
   end type

   type :: BasisCenter
         real(dp) :: coord(3)=(/ZERO,ZERO,ZERO/)
   end type

   type, extends(BasisCenter) :: SolHarmBasisCenter
         real(dp) :: rscale=ONE
         integer :: lmin=0, lmax=-1
         type(SpVecTuple), allocatable :: coeff(:)
   end type

   ! This base class (Basis) should also include a method that yields
   ! the number of unknown expansion coefficients of that basis.
   type :: Basis
         character(80) :: name="unknown"
      contains
         procedure :: get_size => Basis_get_size
   end type

   type, extends(Basis) :: SolHarmBasis
         type(SolHarmBasisCenter), allocatable :: centers(:)
         integer :: solharm_type=MPE_CONST%BASIS_UNDEF
   end type

   type :: DielectricContinuum
         real(dp) :: eps=ONE
         class(Basis), allocatable :: basis
   end type

   type :: InterfacePoint
         real(dp) :: coord(3)
         real(dp) :: normal(3)
         real(dp) :: tangents(3,2) !TODO: currently gets lost upon import
         real(dp) :: area
         real(dp) :: volume
   end type

   type :: DielectricInterface
         ! sampling points (including coordinate system)
         type(InterfacePoint), allocatable :: p(:)
         ! number of boundary conditions at this interface
         integer :: n_bc
         ! indices of dielectric continua at the interface,
         ! dc_ind_pos belongs to the medium in positive normal direction
         integer :: dc_ind_pos, dc_ind_neg
         real(dp) :: surface_energy, volume_energy
   end type

   type :: SpVec
         type(SpVecTuple), allocatable :: a(:)
   end type SpVec

   ! VARIABLES
   ! Note: module variables are initialized in set_mpe_defaults()
   logical, public :: dummy_mpi, module_initialized

   ! data types for MPI communication
   integer :: SpVecTuple_mpi_type
   integer :: InterfacePoint_mpi_type



contains



!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/initialize_mpe_types
!  NAME
!    initialize_mpe_types
!  SYNOPSIS

subroutine initialize_mpe_types(mpi_present)

!  PURPOSE
!    This subroutine creates MPI types for the derived types where necessary.
!  USES
   implicit none
!  ARGUMENTS
   logical, intent(in) :: mpi_present
!  INPUTS
!   o mpi_present -- MPI is present and not just dummy routines
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   if (module_initialized) return

   module_initialized = .true.
   dummy_mpi = .not. mpi_present

   ! create derived MPI data types
   if (.not.dummy_mpi) then
      call SpVecTuple_build_mpi_type(SpVecTuple_mpi_type)
      call InterfacePoint_build_mpi_type(InterfacePoint_mpi_type)
   endif

end subroutine initialize_mpe_types


!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/cleanup_mpe_types
!  NAME
!    cleanup_mpe_types
!  SYNOPSIS

subroutine cleanup_mpe_types()

!  PURPOSE
!    This subroutine frees all MPI types created during initialization.
!  USES
   implicit none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: ierr

   if (.not.module_initialized) return

   module_initialized = .false.

   ! clear derived MPI data types
   if (.not.dummy_mpi) then
      call MPI_TYPE_FREE(SpVecTuple_mpi_type, ierr)
   endif

end subroutine cleanup_mpe_types



!-------------------------------------------------------------------------------
!                                    MPI

!******
!----------------------------------------------------------------------
!****s* mpe_reaction_field/resized_mpi_type
!  NAME
!    resized_mpi_type
!  SYNOPSIS
integer function resized_mpi_type( nblk, blklens, typelist, &
                                   displs, lb, extent )
!  USES
   implicit none
!  PURPOSE
!    Build, commit and return MPI type.
!  ARGUMENTS
   integer, intent(in) :: nblk
   integer, intent(in) :: blklens(nblk), typelist(nblk)
   integer(kind=MPI_ADDRESS_KIND), intent(in) :: displs(nblk), lb, extent
!  INPUTS
!    o n_blocks -- number of blocks
!    o blklens -- size of each block
!    o displacements -- starting address of each block
!    o typelist -- list of existing MPI types
!    o lb -- starting address of derived type
!    o extent -- "real" size of derived type
!  OUTPUT
!    o resized_mpi_type -- new MPI data type
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: ierr, tmp_data_type

   call MPI_Type_create_struct(nblk, blklens, displs-lb, &
                                 typelist, tmp_data_type, ierr)
   call MPI_TYPE_COMMIT(tmp_data_type, ierr)
   call MPI_Type_create_resized(tmp_data_type, lb-lb, extent-lb, &
                                 resized_mpi_type, ierr)
   call MPI_TYPE_FREE(tmp_data_type, ierr)
   call MPI_TYPE_COMMIT(resized_mpi_type, ierr)

end function resized_mpi_type


!******
!----------------------------------------------------------------------
!****s* mpe_reaction_field/InterfacePoint_build_mpi_type
!  NAME
!    InterfacePoint_build_mpi_type
!  SYNOPSIS
subroutine InterfacePoint_build_mpi_type(derived_mpi_type)
!  USES
   implicit none
!  PURPOSE
!    Build MPI type for derived type InterfacePoint which is needed
!    to communicate data of this type via MPI.
!
!    Please note that this treatment requires the items in the derived
!    type to be stored sequentially, so better use the "sequential" keyword
!
!    Customized section consists of blocks of the following type:
!
!    index of block in the derived type
!    iblk = {current number}
!
!    array length of the variable in this block
!    blklens(iblk) = {number of repetitions}
!
!    existing MPI data type for variable type
!    typelist(iblk) = {MPI data type}
!
!    determine starting memory address (displacement) of block
!    call MPI_Get_address(template%{item}, displs(iblk), ierr)
!
!  ARGUMENTS

   integer, intent(out) :: derived_mpi_type

!  INPUTS
!    none
!  OUTPUT
!    o derived_mpi_type -- new MPI data type
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   type(InterfacePoint) :: template(2)
   integer, parameter :: nblk = 4
   integer :: ierr, iblk, blklens(nblk), typelist(nblk)
   integer(kind=MPI_ADDRESS_KIND) :: displs(nblk), lb, extent

   ! CUSTOMIZED

   ! variable: coord
   iblk = 1
   blklens(iblk) = 3
   typelist(iblk) = type_mpi % real_dp
   call MPI_Get_address(template(1)%coord, displs(iblk), ierr)

   ! variable: normal
   iblk = 2
   blklens(iblk) = 3
   typelist(iblk) = type_mpi % real_dp
   call MPI_Get_address(template(1)%normal, displs(iblk), ierr)

   ! variable: tangents
   iblk = 3
   blklens(iblk) = 6
   typelist(iblk) = type_mpi % real_dp
   call MPI_Get_address(template(1)%tangents, displs(iblk), ierr)

   ! variable: area
   iblk = 4
   blklens(iblk) = 1
   typelist(iblk) = type_mpi % real_dp
   call MPI_Get_address(template(1)%area, displs(iblk), ierr)


   ! INVARIANT

   ! get starting address and extent
   call MPI_Get_address(template(1), lb, ierr)
   call MPI_Get_address(template(2), extent, ierr)
   ! use custom function to build and commit type
   derived_mpi_type = resized_mpi_type( nblk, blklens, typelist, &
                                        displs, lb, extent )

end subroutine InterfacePoint_build_mpi_type


!******
!----------------------------------------------------------------------
!****s* mpe_reaction_field/InterfacePoint_vector_mpi_bcast
!  NAME
!    InterfacePoint_vector_mpi_bcast
!  SYNOPSIS
subroutine InterfacePoint_vector_mpi_bcast(v, from_here, mpi_comm)
!  USES
   implicit none
!  PURPOSE
!    Broadcast an array of type InterfacePoint with MPI_Bcast
!    within the MPI communicator.
!  ARGUMENTS
   type(InterfacePoint), intent(inout) :: v(:)
   integer, intent(in) :: from_here, mpi_comm

!  INPUTS
!    o v -- the array to be broadcasted
!    o from_here -- the index of the broadcasting task
!    o mpi_comm -- MPI communicator
!  OUTPUT
!    o v -- overwritten on receiving tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: length, ierr

   if (dummy_mpi) return

   length = size(v)
   call MPI_Bcast(v, length, InterfacePoint_mpi_type, &
                  from_here, mpi_comm, ierr)

end subroutine InterfacePoint_vector_mpi_bcast


!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/SpVecTuple_build_mpi_type
!  NAME
!    SpVecTuple_build_mpi_type
!  SYNOPSIS
subroutine SpVecTuple_build_mpi_type(derived_mpi_type)
!  USES
   implicit none
!  PURPOSE
!    Build MPI type for derived type SpVecTuple which is needed
!    to communicate data of this type via MPI.
!
!    Please note that this treatment requires the items in the derived
!    type to be stored sequentially, so better use the "sequential" keyword
!
!    Customized section consists of blocks of the following type:
!
!    index of block in the derived type
!    iblk = {current number}
!
!    array length of the variable in this block
!    blklens(iblk) = {number of repetitions}
!
!    existing MPI data type for variable type
!    typelist(iblk) = {MPI data type}
!
!    determine starting memory address (displacement) of block
!    call MPI_Get_address(template%{item}, displs(iblk), ierr)
!
!  ARGUMENTS

   integer, intent(out) :: derived_mpi_type

!  INPUTS
!    none
!  OUTPUT
!    o derived_mpi_type -- new MPI data type
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   type(SpVecTuple) :: template(2)
   integer, parameter :: nblk = 2
   integer :: ierr, iblk, blklens(nblk), typelist(nblk)
   integer(kind=MPI_ADDRESS_KIND) :: displs(nblk), lb, extent

   ! CUSTOMIZED

   ! variable: val
   iblk = 1
   blklens(iblk) = 1
   typelist(iblk) = type_mpi % real_dp
   call MPI_Get_address(template(1)%val, displs(iblk), ierr)

   ! variable: ind
   iblk = 2
   blklens(iblk) = 1
   typelist(iblk) = MPI_INTEGER
   call MPI_Get_address(template(1)%ind, displs(iblk), ierr)


   ! INVARIANT

   ! get starting address and extent
   call MPI_Get_address(template(1), lb, ierr)
   call MPI_Get_address(template(2), extent, ierr)
   ! use custom function to build and commit type
   derived_mpi_type = resized_mpi_type( nblk, blklens, typelist, &
                                        displs, lb, extent )

end subroutine SpVecTuple_build_mpi_type


!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/SpVecTuple_vector_mpi_bcast
!  NAME
!    SpVecTuple_vector_mpi_bcast
!  SYNOPSIS
subroutine SpVecTuple_vector_mpi_bcast(sparse, from_here, mpi_comm)
!  USES
   implicit none
!  PURPOSE
!    Broadcast a sparse vector of type SpVecTuple from a given task to all
!    tasks in the communicator.
!  ARGUMENTS

   type(SpVecTuple), allocatable, intent(inout) :: sparse(:)
   integer, intent(in) :: from_here, mpi_comm

!  INPUTS
!    o sparse -- the sparse matrix to be broadcast
!    o from_here -- the index of the broadcasting task
!    o mpi_comm -- MPI communicator
!  OUTPUT
!    o sparse -- set on the receiving tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: nnz, ierr

   if (dummy_mpi) return

   ! first, communicate size of dynamically allocatable arrays
   nnz = 0
   if (allocated(sparse)) &
      nnz = size(sparse)

   call MPI_Bcast(nnz, 1, MPI_INTEGER, &
                  from_here, mpi_comm, ierr)

   ! then, ensure proper allocation if necessary
   if (size(sparse).ne.nnz) then
      if (allocated(sparse)) &
         deallocate(sparse)
      allocate(sparse(nnz))
   endif

   ! finally, broadcast values
   call MPI_Bcast(sparse, nnz, SpVecTuple_mpi_type, &
                  from_here, mpi_comm, ierr)

end subroutine SpVecTuple_vector_mpi_bcast



!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/SpVec_mpi_bcast
!  NAME
!    SpVec_mpi_bcast
!  SYNOPSIS
subroutine SpVec_mpi_bcast(sparse, from_here, mpi_comm)
!  USES
   implicit none
!  PURPOSE
!    Broadcast a sparse vector of type SpVec from a given task to all
!    tasks in the global communicator.
!  ARGUMENTS

   type(SpVec), intent(inout) :: sparse
   integer, intent(in) :: from_here, mpi_comm

!  INPUTS
!    o sparse -- the sparse matrix to be broadcast
!    o from_here -- the index of the broadcasting task
!    o mpi_comm -- MPI communicator
!  OUTPUT
!    o sparse -- set on the receiving tasks
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: nnz, ierr

   if (dummy_mpi) return

   ! first, communicate size of dynamically allocatable arrays
   nnz = 0
   if (allocated(sparse%a)) &
      nnz = size(sparse%a)

   call MPI_Bcast(nnz, 1, MPI_INTEGER, &
                  from_here, mpi_comm, ierr)

   ! then, ensure proper allocation if necessary
   if (size(sparse%a).ne.nnz) then
      if (allocated(sparse%a)) &
         deallocate(sparse%a)
      allocate(sparse%a(nnz))
   endif

   ! finally, broadcast values
   call MPI_Bcast(sparse%a, nnz, SpVecTuple_mpi_type, &
                  from_here, mpi_comm, ierr)

end subroutine SpVec_mpi_bcast



!-------------------------------------------------------------------------------
!              type-bound procedures

!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/Basis_get_size
!  NAME
!    Basis_get_size
!  SYNOPSIS

elemental function Basis_get_size(b) result(s)

!  PURPOSE
!    This function returns the number of basis functions of a basis
!    of type Basis
!
!  USES
   implicit none

!  ARGUMENTS
   class(Basis), intent(in) :: b
   integer :: s
!  INPUTS
!   o b -- basis
!  OUTPUT
!   o s -- number of basis functions
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   s = -1
   select type(b)
   class is (SolHarmBasis)
      s = sum(SolHarmBasisCenter_size_lm(b%centers))
   class default
      ! continue
   end select

end function Basis_get_size


!-------------------------------------------------------------------------------
!              general routines only depending on attributes of types

!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/SolHarmBasisCenter_size_lm
!  NAME
!    SolHarmBasisCenter_size_lm
!  SYNOPSIS

elemental function SolHarmBasisCenter_size_lm(center) result(n)

!  PURPOSE
!    This function returns the number of basis functions of a center
!    of type SolHarmBasisCenter, depending on the lmin and lmax attributes.
!
!  USES
   implicit none

!  ARGUMENTS
   type(SolHarmBasisCenter), intent(in) :: center
   integer :: n
!  INPUTS
!   o center -- basis center
!  OUTPUT
!   o n -- number of basis functions
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   n = (center%lmax+1)**2 - (center%lmin)**2

end function SolHarmBasisCenter_size_lm


!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/SolHarmBasisCenter_lbound_lm
!  NAME
!    SolHarmBasisCenter_lbound_lm
!  SYNOPSIS

elemental function SolHarmBasisCenter_lbound_lm(center) result(lbd)

!  PURPOSE
!    This function returns the last combined (lm) functions index
!    of a center of type SolHarmBasisCenter (depending on the lmin attribute).
!
!  USES
   implicit none

!  ARGUMENTS
   type(SolHarmBasisCenter), intent(in) :: center
   integer :: lbd
!  INPUTS
!   o center -- basis center
!  OUTPUT
!   o ubd -- lower bound of basis functions
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   lbd = (center%lmin)**2 + 1

end function SolHarmBasisCenter_lbound_lm


!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/SolHarmBasisCenter_ubound_lm
!  NAME
!    SolHarmBasisCenter_ubound_lm
!  SYNOPSIS

elemental function SolHarmBasisCenter_ubound_lm(center) result(ubd)

!  PURPOSE
!    This function returns the first combined (lm) functions index
!    of a center of type SolHarmBasisCenter (depending on the lmax attribute).
!
!  USES
   implicit none

!  ARGUMENTS
   type(SolHarmBasisCenter), intent(in) :: center
   integer :: ubd
!  INPUTS
!   o center -- basis center
!  OUTPUT
!   o ubd -- lower bound of basis functions
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   ubd = (center%lmax+1)**2

end function SolHarmBasisCenter_ubound_lm


!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/SpVecTuple_vector_from_dense_vector
!  NAME
!    SpVecTuple_vector_from_dense_vector
!  SYNOPSIS
subroutine SpVecTuple_vector_from_dense_vector( dense, indoffset, threshold, &
                                                sparse )
!  USES
   implicit none
!  PURPOSE
!    Create a sparse vector of type SpVec from a dense vector
!  ARGUMENTS

   real(dp), intent(in) :: dense(:)
   integer, intent(in) :: indoffset
   real(dp), intent(in) :: threshold

   type(SpVecTuple), allocatable, intent(out) :: sparse(:)

!  INPUTS
!    o dense -- dense vector
!    o indoffset -- index offset of first item in dense (0: no offset)
!    o threshold -- minimal absolute value to be considered non-zero
!  OUTPUT
!    o sparse -- vector of sparse tuple
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: i_ind, n_ind, i_sparse, n_sparse

   ! ensure that arrays of sparsevec are not allocated in the beginning
   if (allocated(sparse)) deallocate(sparse)

   n_ind = size(dense)
   n_sparse = 0
   ! count number of non-zero entries
   do i_ind = 1, n_ind
      if (abs(dense(i_ind)).gt.threshold) &
         n_sparse = n_sparse + 1
   enddo ! i_ind

   ! allocate sparse tuple array
   allocate(sparse(n_sparse))

   ! shortcut if possible
   if (n_sparse.eq.0) return

   ! paste values
   i_sparse = 0
   do i_ind = 1, n_ind
      if (abs(dense(i_ind)).gt.threshold) then
         i_sparse = i_sparse + 1
         sparse(i_sparse)%ind = i_ind + indoffset
         sparse(i_sparse)%val = dense(i_ind)
      endif
   enddo ! i_ind

end subroutine SpVecTuple_vector_from_dense_vector


!******
!-------------------------------------------------------------------------------
!****s* mpe_reaction_field/SpVec_from_dense_vector
!  NAME
!    SpVec_from_dense_vector
!  SYNOPSIS
subroutine SpVec_from_dense_vector( dense, indoffset, threshold, sparse )
!  USES
   implicit none
!  PURPOSE
!    Create a sparse vector of type SpVec from a dense vector
!  ARGUMENTS

   real(dp), intent(in) :: dense(:)
   integer, intent(in) :: indoffset
   real(dp), intent(in) :: threshold

   type(SpVec), intent(inout) :: sparse

!  INPUTS
!    o dense -- dense vector
!    o indoffset -- index offset of first item in dense (0: no offset)
!    o threshold -- minimal absolute value to be considered non-zero
!  OUTPUT
!    o sparse -- sparse vector
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: i_ind, n_ind, i_sparse, n_sparse

   ! ensure that arrays of sparsevec are not allocated in the beginning
   if (allocated(sparse%a)) deallocate(sparse%a)

   n_ind = size(dense)
   n_sparse = 0
   ! count number of non-zero entries
   do i_ind = 1, n_ind
      if (abs(dense(i_ind)).gt.threshold) &
         n_sparse = n_sparse + 1
   enddo ! i_ind

   ! allocate sparse tuple array
   allocate(sparse%a(n_sparse))

   ! shortcut if possible
   if (n_sparse.eq.0) return

   ! paste values
   i_sparse = 0
   do i_ind = 1, n_ind
      if (abs(dense(i_ind)).gt.threshold) then
         i_sparse = i_sparse + 1
         sparse%a(i_sparse)%ind = i_ind + indoffset
         sparse%a(i_sparse)%val = dense(i_ind)
      endif
   enddo ! i_ind

end subroutine SpVec_from_dense_vector


!******
!-------------------------------------------------------------------------------
!****s* mpe_types/InterfacePoint_vector_extract
!  NAME
!    InterfacePoint_vector_extract
!  SYNOPSIS
subroutine InterfacePoint_vector_extract(IPs, coords, normals, tangents)
!  USES
   implicit none
!  PURPOSE
!    Create a sparse vector of type SpVec from a dense vector
!  ARGUMENTS
   type(InterfacePoint), intent(in) :: IPs(:)

   real(dp), intent(out), optional :: coords(3,size(IPs))
   real(dp), intent(out), optional :: normals(3,size(IPs))
   real(dp), intent(out), optional :: tangents(3,2,size(IPs))
!  INPUTS
!    o IPs -- vector of InterfacePoint type
!  OUTPUT
!    o coords -- array of coordinates
!    o normals -- array of normal vectors
!    o tangents -- array of tangent vectors
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   integer :: i_p

   if (present(coords)) then
      do i_p = lbound(IPs,1), ubound(IPs,1)
         coords(:,i_p) = IPs(i_p) % coord(:)
      enddo ! i_p
   endif ! coords

   if (present(normals)) then
      do i_p = lbound(IPs,1), ubound(IPs,1)
         normals(:,i_p) = IPs(i_p) % normal(:)
      enddo ! i_p
   endif ! normals

   if (present(tangents)) then
      do i_p = lbound(IPs,1), ubound(IPs,1)
         tangents(:,:,i_p) = IPs(i_p) % tangents(:,:)
      enddo ! i_p
   endif ! tangents

end subroutine InterfacePoint_vector_extract


end module mpe_types

