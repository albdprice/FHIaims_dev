!****h* FHI-aims/mpi_types
!  NAME
!    mpi_types
!  SYNOPSIS

module mpi_types

!  PURPOSE
!    This module provides corresponding MPI types and custom MPI operations 
!    for custom kind types defined in the "types" module.
!
!  USES
   use mpi_headers, only: MPI_DATATYPE_NULL, MPI_OP_NULL
   use types, only: dp

   implicit none

   private
   public :: mpi_types_init, mpi_types_clear

   type :: ttype_mpi
      ! For every type added here, a corresponding initialization needs to
      ! be added in the "init" procedure [Init1]
      ! In order to use reduction operations, you might also want to add
      ! a custom reduction operation [1]
         integer :: real_dp = MPI_DATATYPE_NULL
         integer :: complex_dp = MPI_DATATYPE_NULL
   end type

   type(ttype_mpi), public, protected :: type_mpi


   type :: toperation_mpi
      ! [1] custom operations for the custom MPI data types 
      !  (e.g. for reductions like MPI_Allreduce)
      ! These operations need to be linked to Fortran routines
      ! during the initialization [Init2]
         ! real(dp)
         integer :: sum_real_dp = MPI_OP_NULL
         integer :: max_real_dp = MPI_OP_NULL
         ! complex(dp)
         integer :: sum_complex_dp = MPI_OP_NULL
   end type

   type(toperation_mpi), public, protected :: op_mpi


   logical :: module_initialized = .false.

contains

!------------------------------------------------------------------------------
!
! Custom Operations
!
subroutine mpi_op_sum_real_dp(invec, inoutvec, len, mpitype)
   implicit none
   integer, intent(in) :: len
   integer, intent(in) :: mpitype
   real(dp), intent(in) :: invec(len)
   real(dp), intent(inout) :: inoutvec(len)

   inoutvec = inoutvec + invec
end subroutine

subroutine mpi_op_sum_complex_dp(invec, inoutvec, len, mpitype)
   implicit none
   integer, intent(in) :: len
   integer, intent(in) :: mpitype
   complex(dp), intent(in) :: invec(len)
   complex(dp), intent(inout) :: inoutvec(len)

   inoutvec = inoutvec + invec
end subroutine

subroutine mpi_op_max_real_dp(invec, inoutvec, len, mpitype)
   implicit none
   integer, intent(in) :: len
   integer, intent(in) :: mpitype
   real(dp), intent(in) :: invec(len)
   real(dp), intent(inout) :: inoutvec(len)

   inoutvec = max(inoutvec, invec)
end subroutine


!------------------------------------------------------------------------------
!
! INIT
!
subroutine mpi_types_init()
!  PURPOSE
!    Defines MPI types for custom _basic_ fortran types
!    of a certain <kind>.
!
!    Note:
!     Basic types such as those created by
!        MPI_TYPE_CREATE_F90_INTEGER(R, NEWTYPE, IERROR)
!           [equivalent type to selected_int_kind(R)]
!        MPI_TYPE_CREATE_F90_REAL(P, R, NEWTYPE, IERROR)
!           [equivalent type to selected_real_kind(P,R)]
!        MPI_TYPE_CREATE_F90_COMPLEX(P, R, NEWTYPE, IERROR)
!           [equivalent type to selected_real_kind(P,R)]
!     need not be commited.
!
!  USES
   implicit none

   integer :: ierr

   real(dp) :: var_real_dp
   complex(dp) :: var_complex_dp

   if (module_initialized) return

   ! [Init1] type initializations
   ! real(dp)
   call MPI_TYPE_CREATE_F90_REAL(precision(var_real_dp), range(var_real_dp), &
         type_mpi % real_dp, ierr)
   ! complex(dp)
   call MPI_TYPE_CREATE_F90_COMPLEX(precision(var_complex_dp), range(var_complex_dp), &
         type_mpi % complex_dp, ierr)

   ! [Init2] operations initializations
   ! real(dp)
   call MPI_OP_CREATE(mpi_op_sum_real_dp, .true., op_mpi % sum_real_dp, ierr)
   call MPI_OP_CREATE(mpi_op_max_real_dp, .true., op_mpi % max_real_dp, ierr)
   ! complex(dp)
   call MPI_OP_CREATE(mpi_op_sum_complex_dp, .true., op_mpi % sum_complex_dp, ierr)


   module_initialized = .true.
end subroutine

subroutine mpi_types_clear()
!  PURPOSE
!    Free custom MPI operations and commited types.
!
!  USES
   implicit none

   integer :: ierr

   if (.not.module_initialized) return

   ! [Free1] free commited MPI types
   ! Note: Basic types (such as those created by MPI_TYPE_CREATE_F90_*
   !       need not be commited and thus also not freed

   ! [Free2] free custom MPI operations
   ! real(dp)
   call MPI_OP_FREE(op_mpi % sum_real_dp, ierr)
   call MPI_OP_FREE(op_mpi % max_real_dp, ierr)
   ! complex(dp)
   call MPI_OP_FREE(op_mpi % sum_complex_dp, ierr)


   module_initialized = .false.
end subroutine

end module
