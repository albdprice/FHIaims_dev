

! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_mpi

use mbd_constants, only: dp




implicit none


include 'mpif.h'


private :: dp

interface mpi_all_reduce
    module procedure mpi_all_reduce_real_0d
    module procedure mpi_all_reduce_real_1d
    module procedure mpi_all_reduce_real_2d
end interface

contains

subroutine mpi_all_reduce_real_0d(x, comm)
    real(dp), intent(inout) :: x
    integer, intent(in) :: comm

    real(dp) :: x_buffer
    integer :: ierr

    call MPI_ALLREDUCE(x, x_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    x = x_buffer
end subroutine

subroutine mpi_all_reduce_real(x, n, comm)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n)
    integer, intent(in) :: comm

    real(dp), allocatable :: x_buffer(:)
    integer :: ierr

    allocate (x_buffer(n))
    call MPI_ALLREDUCE(x, x_buffer, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    x = x_buffer
end subroutine

subroutine mpi_all_reduce_real_1d(x, comm)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: comm

    call mpi_all_reduce_real(x, size(x), comm)
end subroutine

subroutine mpi_all_reduce_real_2d(x, comm)
    real(dp), intent(inout) :: x(:, :)
    integer, intent(in) :: comm

    call mpi_all_reduce_real(x, size(x), comm)
end subroutine

integer function mpi_get_rank() result(rank)
    integer :: err

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
end function

end module
