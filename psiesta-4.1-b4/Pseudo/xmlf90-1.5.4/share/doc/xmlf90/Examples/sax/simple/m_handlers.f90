module m_handlers

use xmlf90_sax

private

!
! It defines the routines that are called by the XML parser in response
! to particular events.
!
! In this particular example we just print the names of the elements,
! the attribute list, and the content of the pcdata chunks
!
! A module such as this could use "utility routines" to convert pcdata
! to numerical arrays, and to populate specific data structures.
!
public :: begin_element_handler, end_element_handler, pcdata_chunk_handler

CONTAINS  !=============================================================

subroutine begin_element_handler(name,attributes)
character(len=*), intent(in)   :: name
type(dictionary_t), intent(in) :: attributes

write(unit=*,fmt="(2a)") ">>Begin Element: ", name
write(unit=*,fmt="(a,i2,a)") "--- ", len(attributes), " attributes:"
call print_dict(attributes)
end subroutine begin_element_handler

!--------------------------------------------------
subroutine end_element_handler(name)
character(len=*), intent(in)     :: name

  write(unit=*,fmt="(/,2a)") ">>-------------End Element: ", trim(name)

end subroutine end_element_handler

!--------------------------------------------------
subroutine pcdata_chunk_handler(chunk)
character(len=*), intent(in) :: chunk

write(unit=*,fmt="(a)",advance="no") trim(chunk)

end subroutine pcdata_chunk_handler

end module m_handlers












