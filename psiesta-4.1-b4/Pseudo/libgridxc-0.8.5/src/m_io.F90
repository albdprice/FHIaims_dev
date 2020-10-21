!!@LICENSE
!
!     MODULE m_io
!
! Copyright Alberto Garcia, 1996, 1997, 1998
!
! This module implements an interface to the FORTRAN logical unit
! system. Based on code by Richard Maine.
!
! Alberto Garcia, December 30, 1996
! Rewritten as a single subroutine 
! with multiple entry points, March 7, 1998
! Now hybrid to comply with Siesta "die" interface.
! Converted to a module by J.M.Soler. Aug. 2009
! Convert to F90 by Nick R. Papior, Feb, 2018
!---------------------------------------------------------------
!
MODULE m_io
  !
  !-----------------------------------------------------------------
  !
  !     Used module procedures
  !
  USE sys, only: die   ! Termination routine

  implicit none
  !
  !-----------------------------------------------------------------
  !
  !     Public procedures provided by this module
  !
  PUBLIC :: io_seterr  ! Set standard error unit
  PUBLIC :: io_setout  ! Set standard output unit
  PUBLIC :: io_geterr  ! Get standard error unit
  PUBLIC :: io_getout  ! Get standard output unit
  PUBLIC :: io_assign  ! Get some available IO unit and reserve it
  PUBLIC :: io_reserve ! Reserve a specific IO unit
  PUBLIC :: io_close   ! Close and free a given IO unit
  PUBLIC :: io_status   ! Print all used IO units

  PRIVATE ! Nothing is declared public below this point
  !
  !----------------------------------------------------------------
  !
  !     Module variables
  !
  !     Logical unit management. Units 0 to min_lun-1 are "reserved",
  !     since most of the "typical" files (output, etc) use them.
  !
  !     Logical units min_lun to min_max are managed by this module.
  !
  integer, parameter:: min_lun = 10
  integer, parameter:: max_lun = 99
  integer, parameter:: nunits = max_lun-min_lun+1
  integer, save:: stdout = 6
  integer, save:: stderr = 0
  logical, save:: lun_is_free(min_lun:max_lun) = .true.
  !
  !-----------------------------------------------------------------
  !
  !     Internal and dummy variables
  !
  integer  :: i, iostat
  logical  :: used, named, opened
  character:: filename*50, form*11
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  !
  !     Simple interfaces to modify standard units
  !
  subroutine io_seterr(unit)
    integer,intent(in):: unit
    stderr = unit
  end subroutine io_seterr
  !
  !-----------------------------------------------------------------
  !
  subroutine io_setout(unit)
    integer,intent(in):: unit
    stdout = unit
  end subroutine io_setout
  !
  !-----------------------------------------------------------------
  !
  subroutine io_geterr(unit)
    integer,intent(out):: unit
    unit = stderr
  end subroutine io_geterr
  !
  !-----------------------------------------------------------------
  !
  subroutine io_getout(unit)
    integer,intent(out):: unit
    unit = stdout
  end subroutine io_getout
  !
  !------------------------------------------------------------------     
  !
  !     Logical unit management
  !
  subroutine io_assign(lun)
    integer,intent(out):: lun
    !
    !     Looks for a free unit and assigns it to lun
    !
    do lun= min_lun, max_lun
      if (lun_is_free(lun)) then
        inquire(unit=lun, opened=used, iostat=iostat)
        if (iostat .ne. 0) used = .true.
        lun_is_free(lun) = .false.
        if (.not. used) return
      endif
    enddo
    call die('No luns available in io_assign')

  end subroutine io_assign
  !
  !------------------------------------------------------------------     
  !
  subroutine io_reserve(lun)
    integer,intent(in):: lun
    !
    !     Useful to specify that one needs to use a particular unit number
    !
    !     For example, assume some legacy code expects to work with unit 15:
    !
    !     call io_reserve(15)   ! this call at the beginning of the program
    !     ...
    !     open(15,....)
    !
    inquire(unit=lun, opened=used, iostat=iostat)
    if (iostat .ne. 0) used = .true.
    if (used) call die('Cannot reserve unit. Already connected')
    if (lun .ge. min_lun .and. lun .le. max_lun) &
        lun_is_free(lun) = .false.

  end subroutine io_reserve
  !
  !------------------------------------------------------------------     
  !
  subroutine io_close(lun)
    integer,intent(in):: lun
    !
    !     Use this routine instead of a simple close!!
    !
    close(lun)
    if (lun .ge. min_lun .and. lun .le. max_lun) &
        lun_is_free(lun) = .true.

  end subroutine io_close
  !
  !------------------------------------------------------------------     
  !
  subroutine io_status
    !
    !     Prints a list of the connected logical units and the names of
    !     the associated files
    !

    write(stdout,'(a)') '******** io_status ********'
    do i = 0, max_lun
      inquire(i,opened=opened,named=named,name=filename, &
          form=form,iostat=iostat)
      if (iostat .eq. 0) then
        if (opened) then
          if (named) then
            write(stdout,'(i4,5x,a,5x,a)') i, form, filename
          else
            write(stdout,'(i4,5x,a,5x,a)') i, form, 'No name available'
          endif
        endif
      else
        write(stdout,'(i4,5x,a,5x,a)') i, 'Iostat error'
      endif
    enddo
    write(stdout,'(a)') '********           ********'

  end subroutine io_status

END MODULE m_io

