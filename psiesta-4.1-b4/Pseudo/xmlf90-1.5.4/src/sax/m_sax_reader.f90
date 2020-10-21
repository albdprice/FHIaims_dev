module m_sax_reader

private

integer, parameter, public              :: BUFFER_NOT_CONNECTED = -2048
integer, private, parameter             :: MAXLENGTH = 1024

type, public :: file_buffer_t
private
      logical                           :: connected
      logical                           :: eof
      integer                           :: lun
      character(len=50)                 :: filename
      integer                           :: counter
      character(len=MAXLENGTH)          :: buffer
      integer                           :: line
      integer                           :: col
      integer                           :: pos 
      integer                           :: nchars
      logical                           :: debug
      !
      integer                           :: io_eor
      integer                           :: io_eof
end type file_buffer_t

public  :: get_character, sync_file
public  :: line, column, nchars_processed
public  :: open_file, close_file_buffer, rewind_file, mark_eof_file
public  :: eof_file

private :: fill_buffer
private :: get_unit, find_eor_eof

CONTAINS

!-----------------------------------------
!
subroutine open_file(fname,fb,iostat,record_size,verbose)
character(len=*), intent(in)      :: fname
type(file_buffer_t), intent(out)  :: fb
integer, intent(out)              :: iostat
integer, intent(in), optional     :: record_size
logical, intent(in), optional     :: verbose

iostat = 0

call find_eor_eof(fb%io_eor,fb%io_eof)

fb%connected = .false.

call get_unit(fb%lun,iostat)
if (iostat /= 0) then
   if (fb%debug) print *, "Cannot get unit"
   return
endif

if (present(verbose)) then
   fb%debug = verbose
else
   fb%debug = .false.
endif

if (present(record_size)) then
   open(unit=fb%lun,file=fname,form="formatted",status="old", &
        action="read",position="rewind",recl=record_size,iostat=iostat)
else
   open(unit=fb%lun,file=fname,form="formatted",status="old", &
        action="read",position="rewind",recl=65536,iostat=iostat)
endif
if (iostat /= 0) then
   if (fb%debug) print *, "Cannot open file ", trim(fname), " iostat: ", iostat
   return
endif

fb%connected = .true.
fb%counter = 0
fb%eof = .false.
fb%line = 1
fb%col = 0
fb%filename = fname
fb%pos = 0
fb%nchars = 0
fb%buffer = ""

end subroutine open_file

!-------------------------------------------------
subroutine rewind_file(fb)
type(file_buffer_t), intent(inout)  :: fb

fb%eof = .false.
fb%counter = 0
fb%line = 1
fb%col = 0
fb%pos = 0
fb%nchars = 0
fb%buffer = ""

rewind(unit=fb%lun)

end subroutine rewind_file
!-----------------------------------------
subroutine mark_eof_file(fb)
type(file_buffer_t), intent(inout)  :: fb

fb%eof = .true.

end subroutine mark_eof_file

!-----------------------------------------
subroutine close_file_buffer(fb)
type(file_buffer_t), intent(inout)  :: fb

if (fb%connected) then
    close(unit=fb%lun)
    fb%connected = .false.
endif

end subroutine close_file_buffer

!-------------------------------------------------
function eof_file(fb) result (res)
type(file_buffer_t), intent(in)  :: fb
logical                          :: res

res = fb%eof

end function eof_file
!-----------------------------------------
!-----------------------------------------
! New version, able to cope with arbitrarily long lines 
! (still need to specify a big enough record_size if necessary)
!
subroutine fill_buffer(fb,iostat)
type(file_buffer_t), intent(inout)  :: fb
integer, intent(out)  :: iostat
! 
!
character(len=81)  :: str       ! 80 seems like a good compromise?
                                ! (1 extra for added newline, see below)
integer            :: len=-1    ! initialization for buggy compilers
!
read(unit=fb%lun,iostat=iostat,advance="no",size=len,fmt="(a80)") str

if (iostat == fb%io_eof) then
   
   ! End of file
   if (fb%debug) print *, "End of file."
   return

else if (iostat > 0) then

   ! Hard i/o error
   if (fb%debug) print *, "Hard i/o error. iostat:", iostat
   RETURN

else
!
 if (fb%debug) then
   print *, "Buffer: len, iostat", len, iostat
   print *, trim(str)
 endif

   fb%pos = 0

   if (iostat == 0) then

      !  Normal read, with more stuff left on the line
      !
      fb%buffer = str(1:len) 
      fb%nchars = len

   else         ! (end of record)
      !
      !  End of record. We mark it with an LF, whatever it is the native marker.
      !
!!      fb%buffer = str(1:len) // char(10)
      fb%buffer = str(1:len)             !! Avoid allocation of string
      len = len + 1                      !! by compiler
      fb%buffer(len:len) = char(10)   
      fb%nchars = len
      iostat = 0
   endif

endif

end subroutine fill_buffer

!---------------------------------------------------------------
subroutine get_character(fb,c,iostat)
character(len=1), intent(out) :: c
type(file_buffer_t), intent(inout)  :: fb
integer, intent(out)          :: iostat

character(len=1)   :: c_next

if (.not. fb%connected) then
      iostat = BUFFER_NOT_CONNECTED
      return
endif

if (fb%pos >= fb%nchars) then
      call fill_buffer(fb,iostat)
      if (iostat /= 0) return
endif
fb%pos = fb%pos + 1
c = fb%buffer(fb%pos:fb%pos)
fb%counter = fb%counter + 1              ! Raw counter
fb%col = fb%col + 1
!
! Deal with end-of-line handling on the processor...
!
if (c == char(10)) then
   ! Our own marker for end of line
   fb%line = fb%line + 1
   fb%col = 0
endif
if (c == char(13)) then
   c_next = fb%buffer(fb%pos+1:fb%pos+1)
   if (c_next == char(10)) then
      !
      ! Found CRLF. We replace it by LF, as per specs.
      c = c_next
      fb%pos = fb%pos + 1
      if (fb%debug) print *, "-/-> Removed CR before LF in get_character"
   else
      ! Replace single CR by LF
      c = char(10)
      if (fb%debug) print *, "-/-> Changed CR to LF in get_character -- line++"
      !
   endif
   ! In both cases we increase the line counter and reset the column
   !
   fb%line = fb%line + 1
   fb%col = 0
endif

iostat = 0

end subroutine get_character

!----------------------------------------------------
!----------------------------------------------------
! Error Location functions
!
function line(fb) result (ll)
type(file_buffer_t), intent(in)  :: fb
integer                          :: ll

ll = fb%line
end function line

!----------------------------------------------------
function column(fb) result (col)
type(file_buffer_t), intent(in)  :: fb
integer                          :: col

col = fb%col
end function column
!----------------------------------------------------
!----------------------------------------------------
function nchars_processed(fb) result (nc)
type(file_buffer_t), intent(in)  :: fb
integer                          :: nc

nc = fb%counter
end function nchars_processed
!----------------------------------------------------

subroutine sync_file(fb,iostat)
type(file_buffer_t), intent(inout)  :: fb
integer, intent(out)                :: iostat
!
! Repositions the file so that it matches with
! the stored file_buffer information
!
integer          :: target_counter
character(len=1) :: c

target_counter = fb%counter
call rewind_file(fb)
iostat = 0
do
   if (fb%counter == target_counter) exit
   call get_character(fb,c,iostat)
   if (iostat /= 0) return
enddo

end subroutine sync_file

! ----------------------------------------------------------------------
subroutine get_unit(lun,iostat)

! Get an available Fortran unit number

integer, intent(out)  :: lun
integer, intent(out)  :: iostat

integer :: i
logical :: unit_used

do i = 10, 99
   lun = i
   inquire(unit=lun,opened=unit_used)
   if (.not. unit_used) then
      iostat = 0
      return
   endif
enddo
iostat = -1
lun = -1
end subroutine get_unit
! ----------------------------------------------------------------------

subroutine find_eor_eof(io_eor,io_eof)
!
! Determines the values of the iostat values for End of File and 
! End of Record (in non-advancing I/O)
!
integer, intent(out)           :: io_eor
integer, intent(out)           :: io_eof

integer           :: lun, iostat
character(len=1)  :: c

call get_unit(lun,iostat)

if (iostat /= 0) stop "Out of unit numbers"

open(unit=lun,status="scratch",form="formatted", &
     action="readwrite",position="rewind",iostat=iostat)
if (iostat /= 0)   stop "Cannot open test file"

write(unit=lun,fmt=*)  "a"
write(unit=lun,fmt=*)  "b"

rewind(unit=lun)

io_eor = 0
do
  read(unit=lun,fmt="(a1)",advance="NO",iostat=io_eor) c
  if (io_eor /= 0) exit
enddo

io_eof = 0
do
  read(unit=lun,fmt=*,iostat=io_eof)
  if (io_eof /= 0) exit
enddo

close(unit=lun,status="delete")

end subroutine find_eor_eof

end module m_sax_reader
