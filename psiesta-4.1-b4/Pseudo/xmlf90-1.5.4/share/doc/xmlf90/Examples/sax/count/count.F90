program count
  
  use xmlf90_sax
  use m_count

  integer :: i, iostat
  type(xml_t)  :: fxml

  call open_xmlfile(SAX_DATA_DIR//"big-file.xml",fxml,iostat)
  if (iostat /= 0) stop "Cannot open file."

  call xml_parse(fxml, &
                  begin_element_handler, &
                  end_element_handler,   &
                  pcdata_chunk_handler) 
  
  do i=1,nhash
     write(unit=*,fmt="(3a,i6)") "Number of ",trim(element_hash(i)%elm), & 
                 " elements: ",element_hash(i)%num
  enddo

end program count
