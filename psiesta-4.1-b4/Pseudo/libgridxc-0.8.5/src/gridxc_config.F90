module gridxc_config
!
!  Contains global data for internal operation
!  of the library.
!  This makes it non-thread-safe for now
!
! MPI variables
#ifdef MPI
  integer,             save :: gridxc_comm
  integer,public,      save :: gridxc_myNode= -1, gridxc_totNodes=-1
#else
  integer,public,      save :: gridxc_myNode= 0, gridxc_totNodes=1
#endif
!
! Timing options
!
 logical, public :: gridxc_time_mpi_calls = .false.

CONTAINS

#ifdef MPI
SUBROUTINE gridxc_init(comm, time_mpi_calls)

integer, intent(in)           :: comm   ! NOT optional
! internal operation of the library only
logical, intent(in), optional :: time_mpi_calls

integer :: mpierr
#else
SUBROUTINE gridxc_init()
#endif
  
#ifdef MPI
   gridxc_comm = comm
   call MPI_Comm_Size(comm,gridxc_totNodes,mpierr)
   call MPI_Comm_Rank(comm,gridxc_myNode,mpierr)
   !
   gridxc_time_mpi_calls = .false.
   if (present(time_mpi_calls)) then
      gridxc_time_mpi_calls = time_mpi_calls
   endif

#else
   gridxc_myNode = 0
   gridxc_totNodes = 1
   gridxc_time_mpi_calls = .false.
#endif

end subroutine gridxc_init

end module gridxc_config

!
! This is formally an external routine,
! called only by the instrumented MPI
! routines compiled with the package,
! (and not by other instrumented MPI routines)
!
subroutine gridxc_timer_mpi(prog,iopt)

  use gridxc_config, only: gridxc_time_mpi_calls
  use sys, only: gridxc_timer_start
  use sys, only: gridxc_timer_stop

  character(len=*), intent(in) :: prog
  integer, intent(in)          :: iopt

  if (.not. gridxc_time_mpi_calls) return
  if (iopt == 1) then
     call gridxc_timer_start(prog)
  else if (iopt == 2) then
     call gridxc_timer_stop(prog)
  else
     ! do nothing
  endif
end subroutine gridxc_timer_mpi

