c
      subroutine velect(iter,iconv,id,zelec)
c
      use gridXC, only: gridxc_atomXC, gridxc_setXC
#ifdef LIBXC      
      use xc_f90_types_m
      use xc_f90_lib_m
      use gridXC, only: gridxc_setXC_libXC
#endif      
      use gridXC, only: dp
      use m_libxc_sxc_translation, only: xc_id_t
      use m_libxc_sxc_translation, only: get_xc_id_from_atom_id

      implicit none

c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
      include 'compat.h'
c
c    velect generates the electronic output potential from
c    the electron charge density.  The ionic part is
c    added in dsolv1/dsolv2.
c
c    NOTE:  Velect can be called with zelec=ZEL (all electron)
c                             or with zelec=ZVAL (valence only)
c
C     .. Parameters ..
c
      double precision tiny_charge
      parameter (tiny_charge=1.d-12)
c
      double precision zero, one, pnn
      parameter (zero=0.D0,one=1.D0,pnn=.99D0)
c
C     ..
C     .. Scalar Arguments ..
      double precision zelec
      integer iconv, iter
      character id*1
C     ..
C     .. Arrays in Common ..
c
      double precision s1(nrmax), s2(nrmax), w(3*nrmax), y(nrmax),
     &                 yp(nrmax), ypp(nrmax)
C     ..
C     .. Local Scalars ..
      double precision a1, an, b1, bn, ehart, xlo, xnorm,
     &                 ec, exc, vc, vxc, pi, dx, dc, ex, xccor
      integer i, ierr, isx, ll, nrm, relflag, nspin
#ifdef LIBXC      
      integer           :: nfuncs, x_code, c_code, libxc_family
      integer :: libxc_id(2)
      character(len=120) :: libxc_str
      character(len=60)  :: libxc_type
      external libxc_info
#endif      

      character(len=30) :: xc_type, xc_authors
      integer           :: stat
C     ..
C     .. Local Arrays ..
      double precision dens(nrmax,2), vxcarr(nrmax,2)

C     .. External Subroutines ..
      external ext, splift, spliq, excorr
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      logical leqi
      external leqi
c
C     .. Common blocks ..
      common  s1, s2, w, y, yp, ypp
C     ..

      type(xc_id_t) :: xc_id
c
      pi = 4.d0 * atan(1.d0)
      ehart = zero
c
c      fit cd/r by splines
c
      y(1) = zero
      do 10 i = 2, nr
         y(i) = (cdd(i)+cdu(i))/r(i)
   10 continue
      if (ifcore .eq. 2) then
         do 20 i = 2, nr
            y(i) = y(i) + cdc(i)/r(i)
   20    continue
      end if
      isx = 0
      a1 = zero
      an = zero
      b1 = zero
      bn = zero
      nrm = nr
      call splift(r,y,yp,ypp,nrm,w,ierr,isx,a1,b1,an,bn)
      if (ierr .ne. 1) then
         write(6,9000) ierr
         call ext(420+ierr)
      end if
 9000 format(1x,'****** Error in splift ierr =',i2)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
      xlo = zero
      call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s2,ierr)
      if (ierr .ne. 1) then
         write(6,9010) ierr
         call ext(440+ierr)
      end if
 9010 format(1x,'****** Error in spliq ierr =',i2)
      do 30 i = 1, nr
         ypp(i) = r(i)*ypp(i) + 2*yp(i)
         yp(i) = r(i)*yp(i) + y(i)
         y(i) = r(i)*y(i)
   30 continue
      call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s1,ierr)
      if (ierr .ne. 1) then
         write(6,9020) ierr
         call ext(460+ierr)
      end if
 9020 format(1x,'****** Error in spliq ierr =',i2)
c
c      check normalization
c
      xnorm = zero
      if (zelec .ne. zero) xnorm = zelec/s1(nr)
      if (iter .gt. 3 .and. abs(zelec-s1(nr)) .gt. 0.01D0) then
         if (zelec .lt. s1(nr)+1.0D0) then
            write(6,9030) iter, xnorm
 9030       format(/' warning *** charge density rescaled in',' velect',
     &            /' iteration number',i4,3x,'scaling factor =',f6.3,/)
         else
            xnorm = pnn*xnorm
            write(6,9040) iter, xnorm
 9040       format(/' warning *** charge density partially rescaled in',
     &            ' velect',/' iteration number',i4,3x,
     &            'scaling factor =',f6.3,/)
         end if
      end if
c
c      compute new hartree potential
c      renormalize the charge density
c
      do 40 i = 2, nr
         vod(i) = 2*xnorm*(s1(i)/r(i)+s2(nr)-s2(i))
         vou(i) = vod(i)
         cdd(i) = xnorm*cdd(i)
         cdu(i) = xnorm*cdu(i)
   40 continue
c
c      compute hartree contribution to total energy
c
      if (iconv .eq. 1) then
         ehart = zero
         ll = 4
         do 50 i = 2, nr
            ehart = ehart + ll*(cdd(i)+cdu(i))*vod(i)*rab(i)
            ll = 6 - ll
   50    continue
c
c        Divide by two to account for overcounting of the interactions
c        and by three to finish the Simpson rule calculation...
c
         ehart = 0.5d0*(ehart/3)
      end if
c
c     Add the exchange and correlation potential and calculate
c     the total energy contributions.
c
      if ( leqi(icorr,'gl') .or. leqi(icorr,'hl') .or.
     $     leqi(icorr,'wi') .or. leqi(icorr,'bh') ) then
         use_excorr = .true.
      endif

      if (use_excorr) then

         is_gga = .false.
         call excorr(id,cdd,cdu,cdc,vod,vou,vxc,vc,exc,ec)
         etot(4) = ehart
         etot(5) = vxc
         etot(6) = (3*vc-4*ec)
         etot(7) = exc
c
         return

      else
c
c        New XC scheme based on code from Jose Soler and Carlos Balbas
c
c        Compute dens(i,nspin) = density up, density down
c
         do i=2,nr
            if (ispp .eq. 's') then
               dens(i,1) = cdu(i)/(4.d0*pi*r(i)**2)
               dens(i,2) = cdd(i)/(4.d0*pi*r(i)**2)
            else
               dens(i,1) = 0.5d0*(cdu(i) + cdd(i))/(4.d0*pi*r(i)**2)
               dens(i,2) = dens(i,1)
            endif
            if (ifcore .ge. 1) then
               dens(i,1) = dens(i,1) + 0.5d0 * cdc(i)/(4.d0*pi*r(i)**2)
               dens(i,2) = dens(i,2) + 0.5d0 * cdc(i)/(4.d0*pi*r(i)**2)
            endif
         enddo

c
c        Extrapolate the density at r=0  
c
         dens(1,1) = dens(2,1) - (dens(3,1)-dens(2,1))*r(2)/(r(3)-r(2))
         dens(1,2) = dens(2,2) - (dens(3,2)-dens(2,2))*r(2)/(r(3)-r(2))
         if (dens(1,1) .lt. 0.d0) dens(1,1) = 0.d0
         if (dens(1,2) .lt. 0.d0) dens(1,2) = 0.d0
c
c        Define 'relflag' and 'nspin' for the interface ATOMXC
c
         if (ispp .eq. 'r') relflag = 1
         if (ispp .ne. 'r') relflag = 0
         nspin = 2
c
         r(1) = 0.0d0

         if (icorr == "xc") then
#ifdef LIBXC            
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
           is_gga = .false.
           do i = 1, nfuncs
              call libxc_info(libxc_id(i),libxc_str,libxc_type)
              libxc_family = xc_f90_family_from_id(libxc_id(i))
              if (libxc_family==-1) stop "Bad libxc functional"
              is_gga = is_gga .or. (libxc_family == XC_FAMILY_GGA)
           enddo
           
           call gridxc_setXC_libxc(nfuncs,libxc_id)
#else
           stop "libxc not compiled in"
#endif        
         else
          call get_xc_id_from_atom_id(icorr,xc_id,stat)
          if (stat == 0) then
            xc_type = xc_id%siestaxc_id%family
            xc_authors = xc_id%siestaxc_id%authors
            call gridxc_setxc(1,(/trim(xc_type)/), (/trim(xc_authors)/),
     $                (/1._dp/), (/1._dp/))
          else
            stop 'XC ERROR: invalid functional label'
          endif      
            
          is_gga = (xc_type(1:3) == "GGA" .or. xc_type(1:3) == "VDW")
         endif
       
         call gridxc_atomxc(relflag,nr,nrmax,r,nspin,dens,
     .        ex,ec,dx,dc,vxcarr)      

c
c        Add vxc to total potential and energies
c   
         do i=2,nr
            vou(i) = vou(i) + vxcarr(i,1)
            vod(i) = vod(i) + vxcarr(i,2)
         enddo

clcb
         xccor=0.0d0
         ll = 4
         do i=2,nr
            xccor = xccor + ll * rab(i)*
     .           (vxcarr(i,1)*cdu(i) + vxcarr(i,2)*cdd(i))
            ll = 6 - ll
         enddo

         etot(4) = ehart
         etot(5) = xccor / 3
         etot(6) = xccor - 4*(ex + ec)
         etot(7) = ex + ec 

clcb  
c
c        Obtain total potential at r = 0
c
         vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
         vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))

c     *** lcb-jms modification end ********************
         return

      endif
c
c
      end
#ifdef LIBXC      
!     
!     Adapted from libxc example in wiki
!     Alberto Garcia, for use by PSML module
!
!     Used also by velect.
      
      subroutine libxc_info(id,libxc_str,libxc_type)
      use xc_f90_types_m
      use xc_f90_lib_m

      implicit none

      integer, intent(in) :: id
      character(len=*), intent(out) :: libxc_str, libxc_type

      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
!     integer :: i
      character(len=120) :: s1, s2
!     type(xc_f90_pointer_t) :: str

      call xc_f90_func_init(xc_func, xc_info, id, XC_UNPOLARIZED)

      select case(xc_f90_info_kind(xc_info))
      case(XC_EXCHANGE)
         libxc_type = 'exchange'
      case(XC_CORRELATION)
         libxc_type = 'correlation'
      case(XC_EXCHANGE_CORRELATION)
         write(*, '(a)') 'Correlation'
         libxc_type =  'exchange-correlation'
      case(XC_KINETIC)
         write(*, '(a)') 'Correlation'
         libxc_type =  'kinetic'
      end select

      call xc_f90_info_name(xc_info, s1)
      select case(xc_f90_info_family(xc_info))
      case (XC_FAMILY_LDA);       write(s2,'(a)') "LDA"
      case (XC_FAMILY_GGA);       write(s2,'(a)') "GGA"
      case (XC_FAMILY_HYB_GGA);   write(s2,'(a)') "Hybrid GGA"
      case (XC_FAMILY_MGGA);      write(s2,'(a)') "MGGA"
      case (XC_FAMILY_HYB_MGGA);  write(s2,'(a)') "Hybrid MGGA"
      case (XC_FAMILY_LCA);       write(s2,'(a)') "LCA"
      end select
      write(libxc_str, '(4a)') trim(s1), ' (', trim(s2), ')'

!     i = 0
!     call xc_f90_info_refs(xc_info, i, str, s1)
!     do while(i >= 0)
!     write(*, '(a,i1,2a)') '[', i, '] ', trim(s1)
!     call xc_f90_info_refs(xc_info, i, str, s1)
!     end do

      call xc_f90_func_end(xc_func)

      end subroutine libxc_info
#else

!     Stub   
      subroutine libxc_info(id,libxc_str,libxc_type)

      implicit none

      integer, intent(in) :: id
      character(len=*), intent(out) :: libxc_str, libxc_type

      stop "libxc_info: Libxc not compiled in"
      end subroutine libxc_info
#endif








