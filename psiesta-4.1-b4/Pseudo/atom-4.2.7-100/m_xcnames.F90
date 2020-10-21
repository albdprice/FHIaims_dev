module m_xcnames
!
!  Support for XC functional name translation
!
!  Alberto Garcia, June 10, 2014
!
!  A "-1" code in the SiestaXC field means that the library does not
!  implement the functional 
!
!  A "-2" code in the SiestaXC field means that the functional is
!  available in the old "excorr" library in ATOM.
!
! Actual LibXC codes can be obtained through a call to the function
! xc_f90_functional_get_number in the Fortran interface of libxc. In
! earlier versions of the library, this function was supposed to be
! called without the leading "XC_" in the string...
!
! What we call "libxc_status" here is actually a "status" flag: 
!
!            0: implemented
!           -1: not implemented

private

type, public :: xc_id_t
   character(len=30)  :: siestaxc_type 
   character(len=100) :: siestaxc_authors
   integer            :: siestaxc_code 
   character(len=30)  :: libxc_x 
   character(len=30)  :: libxc_c 
   integer            :: libxc_status
   character(len=2)   :: atom_id
end type xc_id_t


type(xc_id_t), dimension(42) :: xct   ! xct stands for "XC table"
! ^LML

! LDA
data xct(1) / xc_id_t("LDA", "PZ", 1, "XC_LDA_X", "XC_LDA_PZ", 0, "ca") /
data xct(2) / xc_id_t("LDA", "PW92", 2, "XC_LDA_X", "XC_LDA_PW", 0, "pw") /
data xct(3) / xc_id_t("LDA", "Wigner", -2, "XC_LDA_X", "XC_LDA_C_WIGNER", 0, "wi") /
data xct(4) / xc_id_t("LDA", "Hedin-Lundqvist", -2, "XC_LDA_X", "XC_LDA_C_HL", 0, "hl") /
data xct(5) / xc_id_t("LDA", "Gunnarson-Lundqvist", -2, "XC_LDA_X", "XC_LDA_C_GL", 0, "gl") /
data xct(6) / xc_id_t("LDA", "von Barth-Hedin", -2, "XC_LDA_X", "XC_LDA_C_vBH", 0, "bh") /
data xct(7) / xc_id_t("LDA", "CA", 1, "XC_LDA_X", "XC_LDA_PZ", 0, "ca") / ! Alias
data xct(8:10) / 3*xc_id_t("", "", -1, "", "", -1, "") /
! GGA
data xct(11) / xc_id_t("GGA", "PBE", 11, "XC_GGA_X_PBE", "XC_GGA_C_PBE", 0, "pb") /
data xct(12) / xc_id_t("GGA", "B86BPBE", 12, "XC_GGA_X_B86_MGC", "XC_GGA_C_PBE", 0, "bp") /
! ^LML
data xct(13) / xc_id_t("GGA", "B88PBE", 13, "XC_GGA_X_B88", "XC_GGA_C_PBE", 0, "be") /
! ^LML
!"RPBE - Hammer et al"
data xct(14) / xc_id_t("GGA", "RPBE", 14, "XC_GGA_X_RPBE", "XC_GGA_C_PBE", 0, "rp") /
!"revPBE Zhang+Yang"
data xct(15) / xc_id_t("GGA", "revPBE", 15, "XC_GGA_X_PBE_R", "XC_GGA_C_PBE", 0, "rv") /
!"Becke-Lee-Yang-Parr"
data xct(16) / xc_id_t("GGA", "LYP", 16, "XC_GGA_X_B88", "XC_GGA_C_LYP", 0, "bl") / !??
!"Wu-Cohen"
data xct(17) / xc_id_t("GGA", "WC", 17, "XC_GGA_X_WC", "XC_GGA_C_PBE", 0, "wc") / !????
!"Perdew-Burke-Ernzerhof-solid"
data xct(18) / xc_id_t("GGA", "PBEsol", 18, "XC_GGA_X_PBE_SOL", "XC_GGA_C_PBE_SOL", 0, "ps") / 
!
data xct(19) / xc_id_t("GGA", "PBEJsJrLO", 19, "XC_GGA_X_PBE_JSJR", "XC_GGA_C_PBE_???", 0, "jo") / 
data xct(20) / xc_id_t("GGA", "PBEJsJrHEG", 20, "XC_GGA_X_PBE_???", "XC_GGA_C_PBE_???", 0, "jh") / 
data xct(21) / xc_id_t("GGA", "PBEGcGxLO", 21, "XC_GGA_X_PBE_???", "XC_GGA_C_PBE_???", 0, "go") / 
data xct(22) / xc_id_t("GGA", "PBEGcGxHEG", 22, "XC_GGA_X_PBE_???", "XC_GGA_C_PBE_???", 0, "gh") / 
!"Armiento-Mattsson-05"
data xct(23) / xc_id_t("GGA", "AM05", 23, "XC_GGA_X_AM05", "XC_GGA_C_AM05", 0, "am") /
data xct(24) / xc_id_t("GGA", "PW91", 24, "XC_GGA_X_PW91", "XC_GGA_C_PW91", 0, "wp") /
data xct(25:32) / 8*xc_id_t("", "", -1, "", "", -1, "") /
! VDW
!"Dion-et-al--DRSLL"
data xct(33) / xc_id_t("VDW", "DRSLL", 33, "XC_GGA_X_OPTB88_VDW", "XC_VDW_C_DF1", -1, "vw") /
data xct(34) / xc_id_t("VDW", "DRSLL", 34, "XC_GGA_X_OPTB88_VDW", "XC_VDW_C_DF1", -1, "vf") /!alias
!""
data xct(35) / xc_id_t("VDW", "LMKLL", 35, "XC_GGA_X_OPTB88_VDW", "XC_VDW_C_DF2", -1, "vl") /
!""
data xct(36) / xc_id_t("VDW", "KBM", 36, "XC_GGA_X_OPTB88_VDW", "XC_VDW_C_DF1", -1, "vk") /
data xct(37) / xc_id_t("VDW", "C09", 37, "XC_GGA_X_OPTB88_VDW", "XC_VDW_C_DF1", -1, "vc") /
data xct(38) / xc_id_t("VDW", "BH", 38, "XC_GGA_X_CX_VDW", "XC_VDW_C_DF1", -1, "vb") /
data xct(39) / xc_id_t("VDW", "VV", 39, "XC_GGA_X_OPTB88_VDW", "XC_VDW_C_VV10", -1, "vv") /
data xct(40:42) / 3*xc_id_t("", "", -1, "", "", -1, "") /
! ^LML


public :: get_xc_id_from_atom_id, print_xc_id
public :: get_xc_id_from_siestaxc

CONTAINS

  subroutine get_xc_id_from_atom_id(atom_id,xc_id,stat)
    character(len=2), intent(in) :: atom_id
    type(xc_id_t), intent(out)   :: xc_id
    integer, intent(out)         :: stat

    integer :: i

    stat = -1
    do i = 1, size(xct)
       if (xct(i)%atom_id == atom_id) then
          xc_id = xct(i)
          stat = 0
       endif
    enddo
  end subroutine get_xc_id_from_atom_id

  subroutine get_xc_id_from_siestaxc(xc_type,xc_authors,xc_id,stat)
    character(len=*), intent(in) :: xc_type
    character(len=*), intent(in) :: xc_authors
    type(xc_id_t), intent(out)   :: xc_id
    integer, intent(out)         :: stat

    integer :: i

    stat = -1
    do i = 1, size(xct)
       if ((xct(i)%siestaxc_type == xc_type)  .and. &
           (xct(i)%siestaxc_authors == xc_authors))  then
          xc_id = xct(i)
          stat = 0
       endif
    enddo
  end subroutine get_xc_id_from_siestaxc
         
  subroutine print_xc_id(xc_id)
    type(xc_id_t), intent(in)   :: xc_id
    print "(a,'--',a,i3,1x,a,'--',a,i4,1x,a2)",  &
          trim(xc_id%siestaxc_type), &
          trim(xc_id%siestaxc_authors), &
          xc_id%siestaxc_code, &
          trim(xc_id%libxc_x), &
          trim(xc_id%libxc_c), &
          xc_id%libxc_status, &
          trim(xc_id%atom_id)
  end subroutine print_xc_id
    
end module m_xcnames

#ifdef __TEST__
program xcid_test

use m_xcnames

type(xc_id_t) :: xc_id
integer       :: stat
character(len=40) :: id, xc_authors, xc_type

do
 write(*,fmt="(a)",advance="no") "Enter string: "
 read(*,"(a)") id
 if (len_trim(id) == 2) then
    call get_xc_id_from_atom_id(trim(id),xc_id,stat)
    if (stat ==0) then
       call print_xc_id(xc_id)
    else
       print *, "UNKNOWN atom_id"
    endif
 else
    write(*,fmt="(a,a)") "String considered as XC type:", trim(id)
    xc_type = id
    write(*,fmt="(a)",advance="no") "Enter XC authors: "
    read(*,"(a)") xc_authors
    call get_xc_id_from_siestaxc(xc_type,xc_authors,xc_id,stat)
    if (stat ==0) then
       call print_xc_id(xc_id)
    else
       print *, "UNKNOWN SiestaXC id"
    endif
 endif
enddo

end program xcid_test
#endif

 

   
