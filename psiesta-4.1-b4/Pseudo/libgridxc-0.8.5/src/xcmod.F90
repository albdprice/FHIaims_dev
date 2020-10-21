!!@LICENSE
!
!******************************************************************************
! MODULE xcmod
! Stores the information about the XC functional to use,
! and provides routines to set and to get it.
!
!******************************************************************************
! subroutine setXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Sets the xc functional(s) to be used by atomxc and/or cellxc
! ------------------------- INPUT ---------------------------------------------
!     integer,         :: n       ! Number of functionals
!     character(len=*),:: func(n) ! Functional name labels
!     character(len=*),:: auth(n) ! Functional author labels
!     real(dp),        :: wx(n)   ! Functional weights for exchange
!     real(dp),        :: wc(n)   ! Functional weights for correlation
!                                   It should be sum(wx)=sum(wc)=1
!
! Allowed functional/author values:
! XCfunc: 
!   'LDA' or 'LSD' => Local density approximation
!            'GGA' => Generalized gradients approx.
!            'VDW' => Van der Waals functional
! XCauth:
!     'LIBXC-XXXX-OPTIONAL_NAME' => Libxc functional
!
!     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
!           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
!           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
!                     the local density limit of the next:
!            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
!           'RPBE' => GGA Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
!         'revPBE' => GGA Zhang & Yang, PRL 80,890(1998)
!        'B86BPBE' => GGA Becke, JCP 85, 7184 (1986); Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
! ^LML
!        'B88PBE'  => GGA Becke, PRA 38, 3098 (1988); Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
! ^LML
!            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
!             'WC' => GGA Wu-Cohen (see subroutine wcxc)
!         'PBESOL' => GGA Perdew et al, PRL, 100, 136406 (2008)
!           'AM05' => GGA Mattsson & Armiento, PRB, 79, 155101 (2009)
!      'PBEJsJrLO' => GGA Reparametrizations of the PBE functional by
!     'PBEJsJrHEG' => GGA   L.S.Pedroza et al, PRB 79, 201106 (2009) and
!      'PBEGcGxLO' => GGA   M.M.Odashima et al, JCTC 5, 798 (2009)
!     'PBEGcGxHEG' => GGA using 4 different combinations of criteria
!          'DRSLL' => VDW Dion et al, PRL 92, 246401 (2004)
!          'LMKLL' => VDW K.Lee et al, PRB 82, 081101 (2010)
!            'KBM' => VDW optB88-vdW of J.Klimes et al, JPCM 22, 022201 (2010)
!            'C09' => VDW V.R. Cooper, PRB 81, 161104 (2010)
!             'BH' => VDW K. Berland and Per Hyldgaard, PRB 89, 035412 (2014)
!             'VV' => VDW Vydrov-VanVoorhis, JCP 133, 244103 (2010)
!
! ------------------------ USAGE ----------------------------------------------
!   use siestaXC, only: setXC
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!
! --------------------- BEHAVIOUR ---------------------------------------------
! - Stops with an error message if n is larger than internal parameter maxFunc
! - Prints a warning message if sum(wx)/=1 or sum(wc)/=1
!
!******************************************************************************
! subroutine getXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Returns the xc functional(s) that has been previously set
! --------------------- OPTIONAL OUTPUT ---------------------------------------
!     integer         :: n       ! Number of functionals
!     character(len=*):: func(n) ! Functional name labels
!     character(len=*):: auth(n) ! Functional author labels
!     real(dp)        :: wx(n)   ! Functional weights for exchange
!     real(dp)        :: wc(n)   ! Functional weights for correlation
!
! ------------------------ USAGE ----------------------------------------------
!   use precision, only: dp
!   use siestaXC,  only: getXC
!   integer,parameter:: maxFunc = 10
!   character(len=20):: func(maxFunc), auth(maxFunc)
!   real(dp):: wx(maxFunc), wc(maxFunc)
!   call getXC( n, func, auth, wx, wc )
!
! --------------------- BEHAVIOUR ---------------------------------------------
! - Does not change any output array whose size is smaller than nFunc
!
!******************************************************************************

! New, easier to use routines:
! call setxc_family_authors(family,authors)
!
! call setxc_libxc_ids(nfuncs,libxc_ids)
!                   /array of libxc ids/
!
module xcmod

  use precision, only: dp              ! Double precision real kind
  use sys,       only: die             ! Termination routine
  use m_vdwxc,   only: vdw_set_author  ! Sets vdW functional flavour

  implicit none

public:: &
  setXC_family_authors, &! Sets single XC functional in family/author style
  setXC, &! Sets XC functional(s) to be used (comprehensive)
  getXC   ! Returns the XC functional(s) being used

#ifdef LIBXC
public :: setXC_libxc_ids! Sets XC functionals using libxc ids
#endif

private ! Nothing is declared public beyond this point

! These data should be put into a derived type and
! initialized and passed around in a single handle, instead
! of being global

  integer, parameter :: maxFunc = 20
  integer,           save :: nXCfunc=0
  character(len=50), save :: XCauth(MaxFunc), XCfunc(MaxFunc)
  real(dp),          save :: XCweightX(MaxFunc), XCweightC(MaxFunc)

contains

  subroutine setXC_family_authors( family, auth )
    !! Sets a single XC functional in family/author style
    implicit none
    character(len=*),intent(in):: family
    character(len=*),intent(in):: auth

    call setXC(1, [family], [auth], [1.0_dp], [1.0_dp])
  end subroutine setXC_family_authors

  
  subroutine setXC( n, func, auth, wx, wc )
    implicit none
    integer,         intent(in):: n       ! Number of functionals
    character(len=*),intent(in):: func(n) ! Functional name labels
    character(len=*),intent(in):: auth(n) ! Functional author labels
    real(dp),        intent(in):: wx(n)   ! Functl weights for exchng
    real(dp),        intent(in):: wc(n)   ! Functl weights for correl

    integer:: i, j

    if (n>maxFunc) call die('setXC: ERROR: parameter maxFunc too small')
    nXCfunc = n
    XCfunc(1:n) = func(1:n)
    XCauth(1:n) = auth(1:n)
    XCweightX(1:n) = wx(1:n)
    XCweightC(1:n) = wc(1:n)
    do i = 1,n
      if (XCfunc(i)=='VDW' .or. XCfunc(i)=='vdw' .or. XCfunc(i)=='vdW') then
        XCfunc(i) = 'VDW'
        do j = 1,i-1
          if (XCfunc(j)=='VDW' .and. XCauth(j)/=XCauth(i)) &
            call die('setXC ERROR: mixing different VDW authors not allowed')
        end do ! j
        call vdw_set_author( XCauth(i) )
      end if ! (XCfunc(i)=='VDW')
     
      if (XCauth(i)(1:6) == "LIBXC-") then
         call process_libxc_spec(XCfunc(i),XCauth(i))
      endif
      
    end do ! i
   
  end subroutine setXC

  subroutine getXC( n, func, auth, wx, wc )
    implicit none
    integer,         optional,intent(out):: n       ! Number of functionals
    character(len=*),optional,intent(out):: func(:) ! Functional name labels
    character(len=*),optional,intent(out):: auth(:) ! Functional author labels
    real(dp),        optional,intent(out):: wx(:)   ! Functl weights for exchng
    real(dp),        optional,intent(out):: wc(:)   ! Functl weights for correl
    integer:: nf
    nf = nXCfunc
    if (present(n)) n = nf
    if (present(func)) then
       if (size(func)>=nf) func(1:nf) = XCfunc(1:nf)
    end if
    if (present(auth)) then
       if (size(auth)>=nf) auth(1:nf) = XCauth(1:nf)
    end if
    if (present(wx)) then
       if (size(wx)  >=nf) wx(1:nf)   = XCweightX(1:nf)
    end if
    if (present(wc)) then
       if (size(wc)  >=nf) wc(1:nf)   = XCweightC(1:nf)
    end if
  end subroutine getXC

#ifndef LIBXC               
  subroutine process_libxc_spec(func,auth)
    character(len=*), intent(in)    :: func
    character(len=*), intent(inout) ::  auth

    call die("Libxc not compiled in. Cannot handle " //  &
              trim(func) // " " // trim(auth))
  end subroutine process_libxc_spec
  
#else
  
  subroutine process_libxc_spec(func,auth)

    use xc_f90_types_m
    use xc_f90_lib_m

    character(len=*), intent(in)    :: func
    character(len=*), intent(inout) ::  auth

    integer :: iostat, xc_id, idx, xc_id_from_symbol
    character(len=50) :: symbolic_name

         ! Fields are of the form LIBXC-XXXX-SYMBOL
         ! where -SYMBOL is optional if XXXX is a meaningful code
         idx = index(auth(7:),"-")
         if (idx /= 0) then
            ! We have code and symbol fields
            read(auth(7:7+idx-2),iostat=iostat,fmt=*) xc_id
            symbolic_name = auth(7+idx:)
            xc_id_from_symbol = xc_f90_functional_get_number(symbolic_name)
            if (xc_id == 0) then
               ! A zero in the code field signals that we want
               ! to fall back on the symbolic name field 
               if (xc_id_from_symbol < 0) then
                  call die("Cannot get xc_id from " // &
                       trim(symbolic_name))
               else
                  xc_id = xc_id_from_symbol
               endif
            else
               ! Check consistency
               if (xc_id /= xc_id_from_symbol) then
                  call die("Conflicting code field for " // &
                       trim(symbolic_name))
               endif
            endif
            ! Normalize the internal representation 
            write(auth,"(a,i4.4,'-',a)") &
                 "LIBXC-",xc_id, trim(symbolic_name)
         else
            ! Just a code field
            read(auth(7:),iostat=iostat,fmt=*) xc_id
            if (iostat /= 0) call die("Bad libxc code in " &
                                   // trim(auth))
            ! Normalize the internal representation 
            write(auth,"(a,i4.4)") "LIBXC-",xc_id
         endif

         !
         select case (xc_f90_family_from_id (xc_id))
         case (XC_FAMILY_LDA)
            if (func /= "LDA") call die("Family mismatch in " // &
                  trim(func) // " " // trim(auth))
         case (XC_FAMILY_GGA)
            if (func /= "GGA") call die("Family mismatch in " // &
                  trim(func) // " " // trim(auth))
         case default
            call die("Unsupported Libxc family or functional")
         end select
  end subroutine process_libxc_spec

#endif

#ifdef LIBXC
  subroutine setXC_libxc_ids( nfuncs, libxc_ids)
    !! Sets the XC info using libxc numerical codes
    
    use xc_f90_types_m
    use xc_f90_lib_m

    implicit none
    integer, intent(in) :: nfuncs
    !! number of functionals
    integer, intent(in) :: libxc_ids(nfuncs)
    !! numerical libxc codes

    ! automatic arrays
    character(len=10)   :: family(nfuncs), auth(nfuncs)
    real(dp)            :: weight_x(nfuncs), weight_c(nfuncs)
    
    type(xc_f90_pointer_t) :: xc_func, xc_info
    integer :: xc_ispin
    integer :: i

    do i = 1, nfuncs
       !
       ! Determine the kind of functional to assign weights correctly
       !
       xc_ispin = XC_UNPOLARIZED 
       ! 'unpolarized' is the least stringent option: for non-polarized
       ! functionals, the other option might result in an error
       call xc_f90_func_init(xc_func, xc_info, libxc_ids(i), xc_ispin)
       
       select case (xc_f90_info_kind(xc_info))
       case (XC_CORRELATION)
          weight_x(i) = 0.0_dp
          weight_c(i) = 1.0_dp
       case (XC_EXCHANGE)
          weight_x(i) = 1.0_dp
          weight_c(i) = 0.0_dp
       case (XC_EXCHANGE_CORRELATION)
          weight_x(i) = 1.0_dp
          weight_c(i) = 1.0_dp
       case default
          call die("Functional kind not supported")
       end select
       call xc_f90_func_end(xc_func)

       select case (xc_f90_family_from_id (libxc_ids(i)))
       case (XC_FAMILY_LDA)
          family(i) = "LDA"
       case (XC_FAMILY_GGA)
          family(i) = "GGA"
       ! There is probably a (negative) case for bad id ...
       case (-1)
          call die("Bad libxc functional code")
       case default
          family(i) = "other"
       end select
       write(auth(i),"(a,i4.4)") "LIBXC-", libxc_ids(i)
    end do

    if (sum(weight_x(1:nfuncs)) /= 1.0_dp) then
       call die("Wrong exchange weights")
    endif
    if (sum(weight_c(1:nfuncs)) /= 1.0_dp) then
       call die("Wrong correlation weights")
    endif
    
    call setXC(nfuncs, family, auth, weight_x, weight_c)
    
  end subroutine setXC_libxc_ids
# endif /* LIBXC */
  
end module xcmod
