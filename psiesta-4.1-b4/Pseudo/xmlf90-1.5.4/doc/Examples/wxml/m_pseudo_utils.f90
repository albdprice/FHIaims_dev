
      module m_pseudo_utils
!
!     Main module for ps input and output, based on 
!     a derived type representation closely resembling
!     the Froyen data structures.
!
!
!     The radial coordinate is reparametrized to allow a more
!     precise sampling of the area near the origin.

!      r: 0->...
!      x: 0->...
!      i: 1->...

!           r(x) = grid_scale *  [ exp( grid_step*x) - 1 ]
!
!     **** WARNING *****
!     In SIESTA, grid_scale = b, grid_step = a
!     In ATOM,   grid_scale = a, grid_step = b
!     ******************
!
!     pseudo%nr and pseudo%nrval are identical
!     (for ATOM and SIESTA use)
!
!     Working precision should be double precision
!     for backwards binary compatibility 
!
      private

      public :: pseudo_read_formatted
      public :: pseudo_header_print
      public :: pseudo_write_xml
      public :: pseudo_complete
      private :: get_unit

      integer, private, parameter :: dp = selected_real_kind(14)

      type, public  ::  pseudopotential_t
        character(len=2)        :: name
        integer                 :: nr
        integer                 :: nrval
        real(dp)                :: zval
        logical                 :: relativistic
        character(len=2)        :: icorr
        character(len=3)        :: irel
        character(len=4)        :: nicore
        real(dp)                :: grid_scale
        real(dp)                :: grid_step
        character(len=10)       :: method(6)
        character(len=70)       :: text
        integer                 :: npotu
        integer                 :: npotd
        real(dp), pointer       :: r(:)
        real(dp), pointer       :: chcore(:)
        real(dp), pointer       :: chval(:)
        real(dp), pointer       :: vdown(:,:)
        real(dp), pointer       :: vup(:,:)
        integer, pointer        :: ldown(:)
        integer, pointer        :: lup(:)
!
!       Extra fields for more functionality
!
        character(len=10)       :: creator
        character(len=10)       :: date
        character(len=40)       :: flavor
        integer                 :: lmax
        integer, pointer        :: principal_n(:)
        real(dp), pointer       :: occupation(:)
        real(dp), pointer       :: cutoff(:)
      end type pseudopotential_t
!
!     These determine the format for ASCII files
!
      character(len=*), parameter, private ::  &
               fmt_int = "(tr1,i2)" ,                  &
               fmt_nam = "(tr1,a2,tr1,a2,tr1,a3,tr1,a4)", &
               fmt_met = "(tr1,6a10,/,tr1,a70)" ,       &
               fmt_pot= "(tr1,2i3,i5,3es20.12)" ,      &
               fmt_rad = "(4(es20.12))"        ,      &
               fmt_txt = "(tr1,a)"

        CONTAINS

!----
        subroutine pseudo_read_formatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer  :: io_ps, i, j, status
        character(len=70) :: dummy
        real(kind=dp)  :: r2

        call get_unit(io_ps,status)
        if (status /= 0) stop "cannot get unit number"
        open(unit=io_ps,file=fname,form="formatted",status="old",&
             action="read",position="rewind")
        write(6,"(3a)") "Reading pseudopotential information ", &
            "in formatted form from ", trim(fname)


        read(io_ps,fmt=fmt_nam) p%name, p%icorr, p%irel, p%nicore
        read(io_ps,fmt_met) (p%method(i),i=1,6), p%text
        read(io_ps,fmt=fmt_pot) p%npotd, p%npotu, p%nr, &
                                p%grid_scale, p%grid_step, p%zval

        p%nrval = p%nr + 1
        p%nr = p%nr + 1
        allocate(p%r(1:p%nrval))
        read(io_ps,fmt=fmt_txt) dummy
        read(io_ps,fmt=fmt_rad) (p%r(j),j=2,p%nrval)
        p%r(1) = 0.d0
        r2=p%r(2)/(p%r(3)-p%r(2))

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif
        do i=1,p%npotd
           read(io_ps,fmt=fmt_txt) dummy 
           read(io_ps,fmt=fmt_int) p%ldown(i)
           read(io_ps,fmt=fmt_rad) (p%vdown(i,j), j=2,p%nrval)
           p%vdown(i,1) = p%vdown(i,2) - r2*(p%vdown(i,3)-p%vdown(i,2))
        enddo

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
        do i=1,p%npotu
           read(io_ps,fmt=fmt_txt) dummy 
           read(io_ps,fmt=fmt_int) p%lup(i)
           read(io_ps,fmt=fmt_rad) (p%vup(i,j), j=2,p%nrval)
           p%vup(i,1) = p%vup(i,2) - r2*(p%vup(i,3)-p%vup(i,2))
        enddo

        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        read(io_ps,fmt=fmt_txt) dummy
        read(io_ps,fmt=fmt_rad) (p%chcore(j),j=2,p%nrval)
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))

        read(io_ps,fmt=fmt_txt) dummy
        read(io_ps,fmt=fmt_rad) (p%chval(j),j=2,p%nrval)
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))

        close(io_ps)
        end subroutine pseudo_read_formatted
!------

        subroutine vps_init(p)
        type(pseudopotential_t)  :: p
        nullify(p%lup,p%ldown,p%r,p%chcore,p%chval,p%vdown,p%vup)
        end subroutine vps_init

!-------
        subroutine pseudo_header_print(lun,p)
        integer, intent(in) :: lun
        type(pseudopotential_t)  :: p

        integer :: i

        write(lun,"(a)") "<pseudopotential_header>"
        write(lun,fmt=fmt_nam) p%name, p%icorr, p%irel, p%nicore
        write(lun,fmt_met) (p%method(i),i=1,6), p%text
        write(lun,"(a)") "</pseudopotential_header>"

        end subroutine pseudo_header_print
!--------
subroutine pseudo_write_xml(fname,p)
use xmlf90_wxml

character(len=*), intent(in) :: fname
type(pseudopotential_t)      :: p

integer  :: i
type(xmlf_t) :: xf

call xml_OpenFile(fname,xf,indent=.true.)
call xml_AddXMLDeclaration(xf)
call xml_NewElement(xf,"pseudo")
call xml_AddAttribute(xf,"version","0.5")
call xml_NewElement(xf,"header")
call xml_AddAttribute(xf,"symbol",trim(p%name))
!call xml_AddAttribute(xf,"zval",trim(str(p%zval)))
call xml_AddAttribute(xf,"zval",p%zval) ! Overloaded
call xml_AddAttribute(xf,"creator",trim(p%creator))
call xml_AddAttribute(xf,"date",trim(p%date))
call xml_AddAttribute(xf,"flavor",trim(p%flavor))
call xml_AddAttribute(xf,"correlation",trim(p%icorr))

      select case (trim(p%irel))
         case ("isp") 
            call xml_AddAttribute(xf,"relativistic","no")
            call xml_AddAttribute(xf,"polarized","yes")
         case ("nrl")
            call xml_AddAttribute(xf,"relativistic","no")
            call xml_AddAttribute(xf,"polarized","no")
         case ("rel")
            call xml_AddAttribute(xf,"relativistic","yes")
            call xml_AddAttribute(xf,"polarized","no")
      end select
      call xml_AddAttribute(xf,"core-corrections",trim(p%nicore))
call xml_EndElement(xf,"header")

call xml_NewElement(xf,"grid")
      call xml_AddAttribute(xf,"type","log")
      call xml_AddAttribute(xf,"units","bohr")
      call xml_AddAttribute(xf,"scale",p%grid_scale)
      call xml_AddAttribute(xf,"step",p%grid_step)
      call xml_AddAttribute(xf,"npts",p%nr-1)
      call xml_EndElement(xf,"grid")

call xml_NewElement(xf,"semilocal")

      call xml_AddAttribute(xf,"units","rydberg")
      call xml_AddAttribute(xf,"format","r*V")
      call xml_AddAttribute(xf,"npots-down",p%npotd)
      call xml_AddAttribute(xf,"npots-up",p%npotu)

      do i=1,p%npotd
         call xml_NewElement(xf,"vps")
         call xml_AddAttribute(xf,"principal-n", &
                p%principal_n(p%ldown(i)))
         call xml_AddAttribute(xf,"l",p%ldown(i))
         call xml_AddAttribute(xf,"cutoff", &
                p%cutoff(p%ldown(i)))
         call xml_AddAttribute(xf,"occupation", &
                p%occupation(p%ldown(i)))
         call xml_AddAttribute(xf,"spin","-1")

         call xml_NewElement(xf,"radfunc")
         call xml_NewElement(xf,"data")
         call xml_AddArray(xf,p%vdown(i,2:p%nr),fmt_rad)
         call xml_EndElement(xf,"data")
         call xml_EndElement(xf,"radfunc")
         call xml_EndElement(xf,"vps")
      enddo

      do i=1,p%npotu
         call xml_NewElement(xf,"vps")
         call xml_AddAttribute(xf,"principal-n", &
                trim(str(p%principal_n(p%lup(i)))))
         call xml_AddAttribute(xf,"l",trim(str(p%lup(i))))
         call xml_AddAttribute(xf,"cutoff", &
                trim(str(p%cutoff(p%lup(i)))))
         call xml_AddAttribute(xf,"occupation", &
                trim(str(p%occupation(p%lup(i)))))
         call xml_AddAttribute(xf,"spin","-1")

         call xml_NewElement(xf,"radfunc")
         call xml_NewElement(xf,"data")
         call xml_AddArray(xf,p%vup(i,2:p%nr),fmt_rad)
         call xml_EndElement(xf,"data")
         call xml_EndElement(xf,"radfunc")
         call xml_EndElement(xf,"vps")
      enddo

      call xml_EndElement(xf,"semilocal")
      
      call xml_NewElement(xf,"valence-charge")
         call xml_NewElement(xf,"radfunc")
         call xml_NewElement(xf,"data")
         call xml_AddArray(xf,p%chval(2:p%nr),fmt_rad)
         call xml_EndElement(xf,"data")
         call xml_EndElement(xf,"radfunc")
         call xml_EndElement(xf,"valence-charge")

      call xml_NewElement(xf,"pseudocore-charge")
         call xml_NewElement(xf,"radfunc")
         call xml_NewElement(xf,"data")
         call xml_AddArray(xf,p%chcore(2:p%nr),fmt_rad)
         call xml_EndElement(xf,"data")
         call xml_EndElement(xf,"radfunc")
         call xml_EndElement(xf,"pseudocore-charge")

   call xml_EndElement(xf,"pseudo")

   call xml_Close(xf)

end subroutine pseudo_write_xml

!
!  Experimental routine to extract information from "text" field.
!  and to set up more rational fields.
!
subroutine pseudo_complete(p)
type(pseudopotential_t), intent(inout) :: p

integer  :: i, lmax, l, itext, n, status
real(dp) :: zup, zdown, ztot, rc_read

p%creator = p%method(1)
p%date = p%method(2)
p%flavor = p%method(3) // p%method(4) // p%method(5) // p%method(6)

lmax = 0
do i = 1, p%npotd
   lmax = max(lmax,p%ldown(i))
enddo
p%lmax = lmax
allocate(p%principal_n(0:lmax))
allocate(p%occupation(0:lmax))
allocate(p%cutoff(0:lmax))
!
! Decode text into useful information. Assumes l's are increasing from 0
!
if (p%irel=="isp") then
      print *, "Polarized........*************"
      print *, "|", p%text, "|"
      do l=0,min(lmax,3)
         itext=l*17
         read(p%text(itext+1:),iostat=status, &
                     fmt="(i1,tr1,f4.2,tr1,f4.2,tr1,f4.2)") &
                     n, zdown, zup, rc_read
         if (status /=0) STOP "fallo text"
         p%principal_n(l) = n
         p%occupation(l) = zdown+zup
         p%cutoff(l) = rc_read
      enddo
else
      do l=0,min(lmax,3)
         itext=l*17
         read(p%text(itext+1:),iostat=status,fmt="(i1,tr1,f5.2,tr4,f5.2)") &
                      n, ztot, rc_read
         if (status /=0) STOP "fallo text"
         p%principal_n(l) = n
         p%occupation(l) = ztot
         p%cutoff(l) = rc_read
      enddo

endif

end subroutine pseudo_complete


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

end module m_pseudo_utils



