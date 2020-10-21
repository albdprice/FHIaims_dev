program pseudoxml
!
! Converts legacy Froyen-style pseudo files to XML
!
use m_pseudo_utils

type(pseudopotential_t)  :: pseudo

call pseudo_read_formatted("PSF",pseudo)
!
! Complete the information in the data structures if possible
!
call pseudo_complete(pseudo)
call pseudo_header_print(6,pseudo)
!
! Output XML
!
call pseudo_write_xml("PSXML",pseudo)

end program pseudoxml



