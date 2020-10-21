program text
!
! Example of XPATH-lite  processing
!
use xmlf90_xpath

type(xml_t) :: fxml

integer  :: status
character(len=100)  :: title

call open_xmlfile(XPATH_DATA_DIR//"Ba.xml",fxml,status)
if (status /=0) then
   print * , "Cannot open file."
   stop
endif

!call enable_debug(sax=.false.)

!
! Search for and print all "title" elements
!
do
      call get_node(fxml,path="//title",pcdata=title,status=status)
      if (status /= 0)  then
         exit
      else
         print *, "Title found: ", trim(title)
      endif
enddo

end program text














