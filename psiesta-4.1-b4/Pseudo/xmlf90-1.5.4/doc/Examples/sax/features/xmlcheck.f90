program xmlcheck
!
! Checks for well-formedness of an XML file
!
use xmlf90_sax

  integer :: iostat
  type(xml_t)  :: fxml

  call open_xmlfile("INP",fxml,iostat)
  if (iostat /= 0) stop "Cannot open file INP."

  call xml_parse(fxml, verbose = .false.)

  print *, "Characters processed: ", xml_char_count(fxml)

end program xmlcheck














