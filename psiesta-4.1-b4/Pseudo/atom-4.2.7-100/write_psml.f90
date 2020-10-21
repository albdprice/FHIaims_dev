

subroutine write_psml( ray, npotd, npotu, zion )

      use xmlf90_wxml
      use m_libxc_sxc_translation, only: xc_id_t
      use m_libxc_sxc_translation, only: get_xc_id_from_atom_id
      use m_libxc_sxc_translation, only: xc_nfuncs_libxc
      use m_uuid
      use m_common   ! Legacy common blocks

  implicit none

  character(len=*), parameter :: PSML_VERSION = "1.1"
  character(len=*), parameter :: PSML_NAMESPACE = "http://esl.cecam.org/PSML/ns/1.1"
  ! Leave this empty and record only atom's version below
  character(len=*), parameter :: PSML_CREATOR = ""

  integer, parameter :: dp = kind(1.d0)
      
  integer, parameter  :: POLY_ORDER_EXTRAPOL = 7  ! For extrapolation at r=0
  real(dp), parameter :: ryd_to_hartree = 0.5_dp

      type(xmlf_t)  :: xf
      type(xc_id_t) :: xc_id

      integer npotd, npotu
      double precision  :: zion

      character*10      :: polattrib, relattrib, coreattrib
      character*10     :: ray(6)
      character*30 xcfuntype, xcfunparam

      integer                         :: ivps, ip, i
      double precision, allocatable   :: chval(:), f(:)

      double precision                :: total_valence_charge
!
      integer :: stat
      character(len=132) :: line, msg
      character(len=36)  :: uuid
      character(len=32000) :: input_buffer

      character(len=120) :: libxc_string
      character(len=60)  :: libxc_type
      integer :: x_code, c_code
      integer :: nfuncs

      external :: libxc_info

      integer  :: nrl
      real(dp) :: rmax, delta
      integer, allocatable  :: isample(:)
      real(dp), allocatable :: r0(:), f0(:)

      allocate(f(1:nr))


! Digest and dump the information about the pseudopotential flavor
      select case(irel) 

        case('isp') 
          polattrib   = 'yes'
          relattrib   = 'no'

        case('rel') 
          polattrib   = 'no'
          relattrib   = 'dirac'

        case('nrl') 
          polattrib   = 'no'
          relattrib   = 'no'

      end select

! Digest and dump the information about the non-linear core corrections
      select case(nicore) 

        case('pcec', 'fcec', 'fche', 'pche') 
          coreattrib  = 'yes'

        case default
          coreattrib  = 'no'

      end select


! Allocate and define the valence charge density
      allocate(chval(1:nr))

      ! Note that we do not renormalize the charge 
      ! to make it integrate to zion
      do ip = 2, nr
        chval(ip) = (cdd(ip)+cdu(ip))
      enddo

! ---------------------------------------------------------------------
                                                                               
      call xml_OpenFile("PSML",xf, indent=.false.)

      call xml_AddXMLDeclaration(xf,"UTF-8")

      call xml_NewElement(xf,"psml")
      call my_add_attribute(xf,"version",PSML_VERSION)
      call my_add_attribute(xf,"energy_unit","hartree")
      call my_add_attribute(xf,"length_unit","bohr")
      call get_uuid(uuid)
      call my_add_attribute(xf,"uuid",uuid)
      call my_add_attribute(xf,"xmlns",PSML_NAMESPACE)


      call xml_NewElement(xf,"provenance")
      call my_add_attribute(xf,"creator",trim(ray(1)) // PSML_CREATOR)
      call my_add_attribute(xf,"date",ray(2))
      call xml_NewElement(xf,"annotation")
       call my_add_attribute(xf,"action","semilocal-pseudopotential-generation")
      call xml_EndElement(xf,"annotation")
      call xml_NewElement(xf,"input-file")
      call my_add_attribute(xf,"name","INP")
!
!     Note that a file already connected to one unit
!     must not be re-opened with another unit...
!     Since INP is still open at this time, we use
!     INP_COPY (generated in atm.f)
!
      open(44,file="INP_COPY",form="formatted",status="old", &
          position="rewind",action="read")
      input_buffer = ""
      do
         read(44,fmt="(a)",iostat=stat) line
         if (stat .ne. 0) exit
         input_buffer = trim(input_buffer) // trim(line) // char(10) 
      enddo
      close(44)
      call xml_AddCDATASection(xf,trim(input_buffer),.true.)
!
      call xml_EndElement(xf,"input-file")
      call xml_EndElement(xf,"provenance")

        call xml_NewElement(xf,"pseudo-atom-spec")
          call my_add_attribute(xf,"atomic-label",nameat)
          call my_add_attribute(xf,"atomic-number",str(znuc))
          call my_add_attribute(xf,"z-pseudo",str(zion))
          call my_add_attribute(xf,"flavor",ray(3)//ray(4))
          call my_add_attribute(xf,"relativity",relattrib)
          call my_add_attribute(xf,"spin-dft",polattrib)
          call my_add_attribute(xf,"core-corrections",coreattrib)
          
           call xml_NewElement(xf,"annotation")
             call my_add_attribute(xf,"pseudo-energy",&
                                   str(ryd_to_hartree*etot(10)))
           call xml_EndElement(xf,"annotation")

          call do_exchange_correlation(icorr)
          ! This is a child element for now, as
          ! implied in the psml paper.
          call do_configuration(total_valence_charge)

        call xml_EndElement(xf,"pseudo-atom-spec")

!        call get_rmax(mmax,rr,nv,uua,rmax)
        delta = 0.005d0
        rmax = 20.0d0  ! For now
        call get_sampled_grid(nr,r,rmax,delta,nrl,isample,r0)
        allocate(f0(nrl))

        call xml_NewElement(xf,"grid")
          call my_add_attribute(xf,"npts",str(nrl))
          call xml_NewElement(xf,"annotation")
           call my_add_attribute(xf,"type","sampled-log-atom")
           call my_add_attribute(xf,"recipe", &
               "r(i:1..N) = a*(exp(b*(i-1))-1); resampled")
           call my_add_attribute(xf,"recipe-cont","a: scale; b: step")
           call my_add_attribute(xf,"scale",str(a))
           call my_add_attribute(xf,"step",str(b))
           call my_add_attribute(xf,"delta",str(delta))
           call my_add_attribute(xf,"rmax",str(rmax))
          call xml_EndElement(xf,"annotation")


          call xml_NewElement(xf,"grid-data")
            call xml_AddArray(xf,r0(1:nrl))
          call xml_EndElement(xf,"grid-data")
        call xml_EndElement(xf,"grid")

  
        call xml_NewElement(xf,"valence-charge")
          call my_add_attribute(xf,"total-charge", &
              str(total_valence_charge))
          call my_add_attribute(xf,"is-unscreening-charge", &
              "yes")
          call my_add_attribute(xf,"rescaled-to-z-pseudo", &
              "no")
          call xml_NewElement(xf,"radfunc")

            call xml_NewElement(xf,"data")
            call remove_r2(chval(:),r(:),f(:))
            call resample(r,f,nr,r0,isample,f0,nrl)
            call xml_AddArray(xf, force_underflow(f0(1:nrl)))
            call xml_EndElement(xf,"data")
          call xml_EndElement(xf,"radfunc")
        call xml_EndElement(xf,"valence-charge")

        if (coreattrib(1:3) .eq. "yes") then

           call xml_NewElement(xf,"pseudocore-charge")
           call my_add_attribute(xf,"matching-radius",str(rc_core))
           call my_add_attribute(xf,"number-of-continuous-derivatives",&
               str(n_of_continuous_derivs))
           call xml_NewElement(xf,"annotation")
           if (n_of_continuous_derivs == 1) then
              call my_add_attribute(xf,"model-core-charge", &
                  "Original Louie-Froyen-Cohen")
           else
              call my_add_attribute(xf,"model-core-charge", &
                  "Three-parameter Martins fit")
           endif
           call xml_EndElement(xf,"annotation")

           call xml_NewElement(xf,"radfunc")

           call xml_NewElement(xf,"data")
            call remove_r2(cdc(:),r(:),f(:))
            call resample(r,f,nr,r0,isample,f0,nrl)
            call xml_AddArray(xf, force_underflow(f0(1:nrl)))
           call xml_EndElement(xf,"data")
           call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"pseudocore-charge")
        endif

! Down pseudopotentials follow
! (Scalar-relativistic set)
! Cannot output reference energies as these are averages
! over l+1/2 and l-1/2        
! There should be an option to output the 'lj' components
        
        if (npotd > 0) then
           call xml_NewElement(xf,"semilocal-potentials")
           if (relattrib=="dirac") then
              call my_add_attribute(xf,"set","scalar_relativistic")
           else 
              if (polattrib=="yes") then
                 call my_add_attribute(xf,"set","spin_average")
              else
                 call my_add_attribute(xf,"set","non_relativistic")
              endif
           endif

      vpsd: do ivps = 1, lmax
           if (indd(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"slps")
             call my_add_attribute(xf,"n",str(no(indd(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"rc",str(rc(ivps)))
             call my_add_attribute(xf,"flavor",ray(3)//ray(4))

             call xml_NewElement(xf,"radfunc")
               call xml_NewElement(xf,"data")
                 call remove_r(viod(ivps,:),r(:),f(:))
                 call resample(r,f,nr,r0,isample,f0,nrl)
                 call xml_AddArray(xf, ryd_to_hartree * f0(1:nrl))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"slps")
         enddo vpsd

         call xml_EndElement(xf,"semilocal-potentials")
      endif
!
! Up pseudopotentials follow
! 

        if (npotu > 0) then

           call xml_NewElement(xf,"semilocal-potentials")
           if (relattrib=="dirac") then
              call my_add_attribute(xf,"set","spin_orbit")
           else 
              if (polattrib=="yes") then
                 call my_add_attribute(xf,"set","spin_difference")
              endif
           endif

         vpsu: do ivps = 1, lmax
           if (indu(ivps) .eq. 0) cycle
           call xml_NewElement(xf,"slps")
             call my_add_attribute(xf,"n",str(no(indu(ivps))))
             call my_add_attribute(xf,"l",il(ivps))
             call my_add_attribute(xf,"rc",str(rc(ivps)))
             call my_add_attribute(xf,"flavor",ray(3)//ray(4))

             call xml_NewElement(xf,"radfunc")

               call xml_NewElement(xf,"data")
                 call remove_r(viou(ivps,:),r(:),f(:))
                 call resample(r,f,nr,r0,isample,f0,nrl)
                 call xml_AddArray(xf, ryd_to_hartree * f0(1:nrl))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"slps")
        enddo vpsu
        call xml_EndElement(xf,"semilocal-potentials")
        endif

! Dump of the pseudowave functions
        call xml_NewElement(xf,"pseudo-wave-functions")
        if (relattrib=="dirac") then
           call my_add_attribute(xf,"set","scalar_relativistic")
        else
           if (polattrib=="yes") then
              call my_add_attribute(xf,"set","spin_average")
           else
              call my_add_attribute(xf,"set","non_relativistic")
           endif
        endif
        ! We only write the "test" wavefunctions ("major" = "average")
        pswfd: do i = 1, nshells_stored
           call xml_NewElement(xf,"pswf")
             call my_add_attribute(xf,"n",str(n_pswf(i)))
             call my_add_attribute(xf,"l",il(l_pswf(i)+1))
             call my_add_attribute(xf,"energy_level", &
                     str(ryd_to_hartree*e_pswf(i)))
             call xml_NewElement(xf,"radfunc")

               call xml_NewElement(xf,"data")
                 call remove_r(pswf(i,:),r(:),f(:))
                 call resample(r,f,nr,r0,isample,f0,nrl)
                 call xml_AddArray(xf, force_underflow(f0(1:nrl)))
               call xml_EndElement(xf,"data")
             call xml_EndElement(xf,"radfunc")
           call xml_EndElement(xf,"pswf")
        enddo pswfd

        call xml_EndElement(xf,"pseudo-wave-functions")


        call xml_EndElement(xf,"psml")
      call xml_Close(xf)

      deallocate(chval)

      CONTAINS

      subroutine do_configuration(total_valence_charge)
      integer, parameter :: dp = selected_real_kind(10,100)

      real(dp), intent(out) :: total_valence_charge
      integer  :: i, lp      
      real(dp) :: occ_down, occ_up, occupation

      call  get_total_valence_charge(total_valence_charge)

      call xml_NewElement(xf,"valence-configuration")

        ! this call is needed here before creating any sub-elements...
        call my_add_attribute(xf,"total-valence-charge", &
                str(total_valence_charge))

      i = ncore
      do 
         i = i + 1
         if (i > norb) exit
         lp = lo(i) + 1
         occ_down = zo(i)
         occ_up   = 0.0_dp
         if (is_split_shell(lp)) then
            i = i + 1
            occ_up = zo(i)
         endif
         occupation = occ_down + occ_up
         if (occupation .lt. 1.0e-10_dp) cycle
         call xml_NewElement(xf,"shell")
         call my_add_attribute(xf,"n",str(no(i)))
         call my_add_attribute(xf,"l",il(lp))
         call my_add_attribute(xf,"occupation",str(occupation))
         if (polarized .and. is_split_shell(lp)) then
            call my_add_attribute(xf,"occupation-down",str(occ_down))
            call my_add_attribute(xf,"occupation-up",str(occ_up))
         endif
         call xml_EndElement(xf,"shell")
      enddo
      call xml_EndElement(xf,"valence-configuration")
      end subroutine do_configuration

      subroutine do_exchange_correlation(icorr)
      character(len=*), intent(in) :: icorr

      integer :: nfuncs, x_code, c_code, libxc_family
      integer :: libxc_id(2)
      character(len=120) :: libxc_str
      character(len=60)  :: libxc_type
      external libxc_info

      ! Digest and dump the information about
      ! the exchange and correlation functional

      call xml_NewElement(xf,"exchange-correlation")

      if (icorr == "xc") then
                    ! libxc functional specification
           ! xc_code is of the form YYYYXXXX, or just ....XXXX
           if (xc_code < 10000) then
            ! Special syntax for single functional
            ! (For example, Teter exch-corr functional: xc_code=0020
            nfuncs = 1
            libxc_id(1) = xc_code
           else
            x_code = xc_code/10000
            c_code = xc_code - 10000*x_code
            nfuncs = 2 
            libxc_id = (/ x_code, c_code /)
           endif

          call xml_NewElement(xf,"annotation")
          call my_add_attribute(xf,"atom-xc-code",icorr)
          call my_add_attribute(xf,"atom-libxc-code",str(xc_code))
          call xml_EndElement(xf,"annotation")

          call xml_NewElement(xf,"libxc-info")
          call my_add_attribute(xf,"number-of-functionals", &
                               str(nfuncs))
          do i = 1, nfuncs
           call xml_NewElement(xf,"functional")
           call libxc_info(libxc_id(i),libxc_string,libxc_type)
           call my_add_attribute(xf,"name",trim(libxc_string))
           call my_add_attribute(xf,"type",trim(libxc_type))
           call my_add_attribute(xf,"id",str(libxc_id(i)))
           call xml_EndElement(xf,"functional")
        enddo
        call xml_EndElement(xf,"libxc-info")

      else
          call get_xc_id_from_atom_id(icorr,xc_id,stat)
          if (stat /= 0) then
             stop "Wrong icorr code!"
          endif
          call xml_NewElement(xf,"annotation")
          call my_add_attribute(xf,"atom-xc-code",icorr)
          call my_add_attribute(xf,"siesta-xc-type", &
                               trim(xc_id%siestaxc_id%family))
          call my_add_attribute(xf,"siesta-xc-authors", &
                               trim(xc_id%siestaxc_id%authors))
          call xml_EndElement(xf,"annotation")

          call xml_NewElement(xf,"libxc-info")
          call my_add_attribute(xf,"number-of-functionals", &
                               str(xc_nfuncs_libxc(xc_id)))
           do i = 1, xc_nfuncs_libxc(xc_id)
              call xml_NewElement(xf,"functional")
               call my_add_attribute(xf,"name", &
                           trim(xc_id%libxc_id(i)%name))
               call my_add_attribute(xf,"type", &
                           trim(xc_id%libxc_id(i)%xc_kind%str))
               call my_add_attribute(xf,"id", &
                           str(xc_id%libxc_id(i)%code))
              call xml_EndElement(xf,"functional")
           enddo
          call xml_EndElement(xf,"libxc-info")
!
       endif
       
       call xml_EndElement(xf,"exchange-correlation")
       
      end subroutine do_exchange_correlation
      
      subroutine get_total_valence_charge(tot_occ)
      integer, parameter :: dp = selected_real_kind(10,100)
      real(dp), intent(out) :: tot_occ 

      integer :: i, lp
      real(dp):: occ_down, occ_up

      tot_occ = 0.0_dp

      i = ncore
      do 
         i = i + 1
         if (i > norb) exit
         lp = lo(i) + 1
         occ_down = zo(i)
         occ_up   = 0.0_dp
         if (is_split_shell(lp)) then
            i = i + 1
            occ_up = zo(i)
         endif
         tot_occ = tot_occ + occ_down + occ_up
      enddo
      end subroutine get_total_valence_charge

      logical function is_split_shell(lp)
      integer, intent(in) :: lp


      is_split_shell = .false.
      if (polarized) then
         is_split_shell = .true.
      else if (relativistic) then
         if (lp /= 1) is_split_shell = .true.
      endif
      end function is_split_shell
         
      subroutine my_add_attribute(xf,name,value)
      type(xmlf_t), intent(inout)   :: xf
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: value

       call xml_AddAttribute(xf,name,trim(value))
      end subroutine my_add_attribute

      subroutine remove_r(rf,r,f)
      ! Removes a factor of r from rf to get f
      ! At the same time, it extrapolates to zero (r(1))

      double precision, intent(in)  :: rf(:), r(:)
      double precision, intent(out) :: f(:)
      
      integer i
      double precision :: r2

      do i = 2, nr
         f(i) = rf(i)/r(i)
      enddo

      r2 = r(2)/(r(3)-r(2))
      f(1) = f(2) - (f(3)-f(2))*r2

      end subroutine remove_r

      subroutine remove_r2(rf,r,f)
      ! Removes a factor of r**2 from rf to get f
      ! At the same time, it extrapolates to zero (r(1))

      double precision, intent(in)  :: rf(:), r(:)
      double precision, intent(out) :: f(:)
      
      integer i
      double precision :: r2

      do i = 2, nr
         f(i) = rf(i)/(r(i)*r(i))
      enddo

      r2 = r(2)/(r(3)-r(2))
      f(1) = f(2) - (f(3)-f(2))*r2

      end subroutine remove_r2

      pure elemental function force_underflow(x)
      double precision, intent(in) ::  x
      double precision  force_underflow

!     Avoid very small numbers that might need a three-character
!     exponent field in formatted output
      
      if (abs(x) .lt. 1.0d-99) then
         force_underflow = 0.0d0
      else
         force_underflow = x
      endif

      end function force_underflow

  subroutine get_sampled_grid(mmax,rr,rmax,delta,nrl,isample,r0)
   integer, intent(in)   :: mmax
   real(dp), intent(in)  :: rr(:)
   real(dp), intent(in)  :: rmax
   real(dp), intent(in)  :: delta
   integer, intent(out)  :: nrl
   integer, allocatable, intent(out) :: isample(:)
   real(dp), allocatable, intent(out) :: r0(:)

   integer  :: is, j
   real(dp) :: rs

   ! First scan to get size of sampled grid
   is = 1
   rs = 0.0_dp
   do j = 1, mmax
      if (rr(j) > rmax) exit
      if ((rr(j)-rs) < delta) cycle
      is = is + 1
      rs = rr(j)
   enddo
   
   nrl = is
   allocate(isample(nrl),r0(nrl))

   is = 1
   r0(is) = 0.0_dp
   isample(is) = 0
   do j = 1, mmax
      if (rr(j) > rmax) exit
      if ((rr(j)-r0(is)) < delta) cycle
      is = is + 1
      r0(is) = rr(j)
      isample(is) = j
   enddo
 end subroutine get_sampled_grid

  subroutine resample(rr,ff,mmax,r0,isample,f0,nrl)
   integer, intent(in)   :: mmax
   real(dp), intent(in)  :: rr(:)
   real(dp), intent(in)  :: ff(:)
   integer, intent(in)   :: nrl
   real(dp), intent(in)  :: r0(:)
   integer, intent(in)   :: isample(:)
   real(dp), intent(out) :: f0(:)

   integer  :: is
   real(dp) :: val
   
   do is = 2, nrl
      f0(is) = ff(isample(is))
   enddo
   ! Choice of treatments of point at r=0
   ! Polynomial extrapolation with sampled points
   call dpnint1(POLY_ORDER_EXTRAPOL,r0(2:),f0(2:),nrl-1,0.0_dp,val,.false.)
   f0(1) = val
   ! Simply set f0(r=0) = ff(r=r1)
   !...
   ! Others
   ! ...
   
 end subroutine resample

!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! 
! Modified by Alberto Garcia, March 2015
! This routine is included in this module with permission from D.R. Hamann.
!
 subroutine dpnint1(npoly, xx, yy, nn, r, val, debug)

! Modified by Alberto Garcia, March 2015 from routine
! dpnint by D.R. Hamann. 
! Changes:
!   -- A single value is returned
!   -- It can extrapolate, instead of stopping,
!      when called with an abscissa outside the
!      data range.
!   -- If the number of data points is less than
!      npoly+1, npoly is implicitly reduced, without
!      error, and without warning.
!   -- Debug interface 
!
! local polynomial interpolation of data yy on nn points xx
! giving value val on point r
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output interpolated value val on point r

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp), intent(in) :: xx(*),yy(*)
 real(dp), intent(in) :: r
 real(dp), intent(out) :: val
 integer, intent(in)   ::  nn,npoly
 logical, intent(in)   ::  debug

!Local variables
 real(dp) :: sum,term,zz
 integer ii,imin,imax,iprod,iy,istart,kk,iend

! interval halving search for xx(ii) points bracketing r

   imin = 1
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(r>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do


   zz=r

!   if (debug) print *, "imin, imax: ", imin, imax

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)
   iend = min(istart+npoly,nn)

 !  if (debug) print *, "istart, iend: ", istart, iend
   sum=0.0d0
   do iy=istart,iend
    if(yy(iy)==0.0d0) cycle
    term=yy(iy)
    do iprod=istart, iend
     if(iprod==iy) cycle
     term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
    end do
    sum=sum+term
   end do
   val=sum

 end subroutine dpnint1


end subroutine write_psml



