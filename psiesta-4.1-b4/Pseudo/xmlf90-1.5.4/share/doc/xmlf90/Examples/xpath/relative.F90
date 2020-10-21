program relative
!
! Example of XPATH-lite  processing
!
use xmlf90_xpath

type(dictionary_t) :: attributes
type(xml_t) :: fxml

integer  :: status

call open_xmlfile(XPATH_DATA_DIR//"Ba.xml",fxml,status)
if (status /=0) then
   print * , "Cannot open file."
   stop
endif

!call enable_debug(sax=.false.)

!
job_search: do
   !
   ! This will search for all the 'job' elements and all the
   ! 'shell' elements with l=0 contained in them at any depth
   ! (relative search).

   call mark_node(fxml,path="/atom/job",attributes=attributes,status=status)
   if (status /= 0)  then
      print *, "No more 'job' elements"
      exit job_search
   else
      print *, ">>>>>>>>>>> New job: "
      call print_dict(attributes)
   endif

   shell_search: do
      !
      ! The initial dot (.) signals a relative search
      !
      call get_node(fxml,path=".//shell",att_name="l", &
           att_value="0",attributes=attributes,status=status)
      if (status /= 0)  then
         print *, "end of job"
         exit shell_search
      endif
      print *, " Found Shell with l=0: "
      call print_dict(attributes)
      print *, "------------------------------------***"
   enddo shell_search

enddo job_search

end program relative














