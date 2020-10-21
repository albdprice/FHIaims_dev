program pseudo_read
!
! Example of XPATH-lite  processing for pseudo xml file
! Shows the use of constrained searches, context delegation, etc.
!
use xmlf90_xpath
use m_pseudo_types

type(dictionary_t) :: attributes
type(xml_t) :: fxml

type(pseudo_t), target, save :: pseudo
type(grid_t),  save          :: global_grid
!
! Pointers to make it easier to manage the data
!
type(header_t),  pointer   :: hp
type(vps_t),  pointer      :: pp

integer               :: status, ndata
character(len=200)    :: value

!-----------------------------------------------------------------
call open_xmlfile(XPATH_DATA_DIR//"pseudo.xml",fxml,status)
if (status /=0) call die("Cannot open file.")

!call enable_debug(sax=.false.)

!
!------------------------------------------------------------
! Root element with version information
!
call get_node(fxml,path="/pseudo",attributes=attributes,status=status)
if (status /= 0)  call die("Cannot find pseudo element")

         call get_value(attributes,"version",value,status)
         if (value == "0.5") then
            print *, "Processing a PSEUDO version 0.5 XML file"
         else
            call die("Can only work with PSEUDO version 0.5 XML files")
         endif

!------------------------------------------------------------
! Header
!
call get_node(fxml,path="/pseudo/header", &
              attributes=attributes,status=status)
if (status /= 0)  call die("Cannot find /pseudo/header")

         hp => pseudo%header
         
         call get_value(attributes,"symbol",hp%symbol,status)
         if (status /= 0 ) call die("Cannot determine atomic symbol")

         call get_value(attributes,"zval",value,status)
         if (status /= 0 ) call die("Cannot determine zval")
         read(unit=value,fmt=*) hp%zval
!
         call get_value(attributes,"creator",hp%creator,status)
         if (status /= 0 ) hp%creator="unknown"

         call get_value(attributes,"flavor",hp%flavor,status)
         if (status /= 0 ) hp%flavor="unknown"

         call get_value(attributes,"relativistic",value,status)
         if (status /= 0 ) value = "no"
         hp%relativistic = (value == "yes")

         call get_value(attributes,"polarized",value,status)
         if (status /= 0 ) value = "no"
         hp%polarized = (value == "yes")

         call get_value(attributes,"core-corrections", &
                                    hp%core_corrections,status)
         if (status /= 0 ) hp%core_corrections = "nc"


!------------------------------------------------------------
! Global grid information
!
call rewind_xmlfile(fxml)
call get_node(fxml,path="/pseudo/grid", &
              attributes=attributes,status=status)

if (status == 0)  then
   print *, "This file has a global grid... "
   call get_grid_data(attributes,global_grid)
else
   global_grid%npts = 0          ! To flag absence of global grid info
endif
!
!------------------------------------------------------------
! Valence charge
!
call rewind_xmlfile(fxml)
!
call mark_node(fxml,path="/pseudo/valence-charge", &
              attributes=attributes,status=status)
if (status == 0)  then
   !
   !  Get the data (and possible private grid)
   !
   call get_radfunc_data(fxml,global_grid,pseudo%valence_charge)
endif
!
!------------------------------------------------------------
! Core charge
!
call rewind_xmlfile(fxml)
!
call mark_node(fxml,path="/pseudo/pseudocore-charge", &
              attributes=attributes,status=status)
if (status == 0)  then
   !
   !  Get the data (and possible private grid)
   !
   call get_radfunc_data(fxml,global_grid,pseudo%core_charge)
endif
!
!------------------------------------------------------------
! Semilocal pseudopotentials
!
call rewind_xmlfile(fxml)
!
call get_node(fxml,path="//semilocal", &
              attributes=attributes,status=status)
if (status /= 0)  call die("Cannot find semilocal element")

         call get_value(attributes,"npots-down",value,status)
         if (status /= 0 ) call die("Cannot determine npots-down")
         read(unit=value,fmt=*) pseudo%npots_down

         call get_value(attributes,"npots-up",value,status)
         if (status /= 0 ) call die("Cannot determine npots-up")
         read(unit=value,fmt=*) pseudo%npots_up

!
! Loop over pseudopotentials
!
pseudo%npots = 0
do
   !
   ! This will search for all the 'vps' elements, marking the context 
   ! in turn
   !
      call mark_node(fxml,path="//vps",attributes=attributes,status=status)
      if (status /= 0) exit          ! exit loop

         pseudo%npots = pseudo%npots + 1
         pp => pseudo%pot(pseudo%npots)

         call get_value(attributes,"l",value,status)
         if (status /= 0 ) call die("Cannot determine l for Vps")
         read(unit=value,fmt=*) pp%l

         call get_value(attributes,"principal-n",value,status)
         if (status /= 0 ) call die("Cannot determine n for Vps")
         read(unit=value,fmt=*) pp%n

         call get_value(attributes,"cutoff",value,status)
         if (status /= 0 ) call die("Cannot determine cutoff for Vps")
         read(unit=value,fmt=*) pp%cutoff

         call get_value(attributes,"occupation",value,status)
         if (status /= 0 ) call die("Cannot determine occupation for Vps")
         read(unit=value,fmt=*) pp%occupation

         call get_value(attributes,"spin",value,status)
         if (status /= 0 ) call die("Cannot determine spin for Vps")
         read(unit=value,fmt=*) pp%spin

         !
         !  Get the data (and possible private grid)
         !
         call get_radfunc_data(fxml,global_grid,pp%V)
         !
         ! After context delegation it is essential to sync the handle
         ! (or to rewind it)
         !
         call sync_xmlfile(fxml,status)
enddo

!
!  Show some of the information
!
call dump_pseudo(pseudo)

!=======================================================================
CONTAINS

!-----------------------------------------------------------------------
subroutine get_radfunc_data(context,global_grid,rp)
!
! Example of routine which packages parsing functionality for a
! common element. The <radfunc> element can appear under <vps>,
! <valence-charge>, and <pseudocore-charge> elements. 
! In all cases the parsing steps are exactly  the same.
! This routine accepts the appropriate context handle and returns
! the data structure.
!
type(xml_t), intent(in)      :: context
type(grid_t), intent(in)     :: global_grid
type(radfunc_t), intent(out) :: rp

type(xml_t)           :: ff
character(len=2000)   :: pcdata

ff = context           ! It inherits the "ancestor element" markings, etc

      call get_node(ff,path="./radfunc/grid", &
              attributes=attributes,status=status)
      if (status == 0)  then
         print *, " >> local grid found"
         call get_grid_data(attributes,rp%grid)
      else
         rp%grid = global_grid
      endif

      ff = context
      call sync_xmlfile(ff,status)  ! Go back to beginning of context 

      call get_node(ff,path="./radfunc/data", &
              pcdata=pcdata,status=status)
      if (status < 0) call die("Cannot find data element")
      if (status > 0) call die("Not enough space for pcdata")
      if (rp%grid%npts == 0) call die("Need grid information!")
      allocate(rp%data(rp%grid%npts))
      ndata = 0             ! To start the build up
      call build_data_array(pcdata,rp%data,ndata)
      if (ndata /= size(rp%data)) STOP "npts mismatch"
end subroutine get_radfunc_data
!-----------------------------------------------------------------------
subroutine get_grid_data(attributes,grid)
type(dictionary_t), intent(in)  :: attributes
type(grid_t), intent(out)       :: grid

         call get_value(attributes,"type",grid%type,status)
         if (status /= 0 ) call die("Cannot determine grid type")

         call get_value(attributes,"npts",value,status)
         if (status /= 0 ) call die("Cannot determine grid npts")
         read(unit=value,fmt=*) grid%npts

         call get_value(attributes,"scale",value,status)
         if (status /= 0 ) call die("Cannot determine grid scale")
         read(unit=value,fmt=*) grid%scale

         call get_value(attributes,"step",value,status)
         if (status /= 0 ) call die("Cannot determine grid step")
         read(unit=value,fmt=*) grid%step

end subroutine get_grid_data

!-----------------------------------------------------------------------
      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) then
         write(unit=0,fmt="(a)") trim(str)
      endif
      write(unit=0,fmt="(a)") "Stopping Program"
      stop
      end subroutine die

end program pseudo_read













