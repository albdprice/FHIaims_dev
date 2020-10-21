program example
!
! Example driver for a stand-alone parsing of an xml document
!
use xmlf90_sax
use m_handlers      ! Defines begin_element, end_element, pcdata_chunk, etc

  integer :: iostat
  type(xml_t)  :: fxml

  call open_xmlfile(SAX_DATA_DIR//"test.xml",fxml,iostat)
  if (iostat /= 0) stop "Cannot open file."

  call xml_parse(fxml, &
               begin_element_handler = begin_element_handler , &
               end_element_handler = end_element_handler, &
               pcdata_chunk_handler = pcdata_chunk_handler, &
               comment_handler = comment_handler, &
               xml_declaration_handler = xml_declaration_handler, &
               sgml_declaration_handler = sgml_declaration_handler, &
               verbose = .false., &
               empty_element_handler = empty_element_handler)

end program example














