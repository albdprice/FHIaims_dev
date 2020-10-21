!!****f* FHI-aims/ext_elsi_yes_no
!!
!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides the function ext_elsi_yes_no, which returns either true or
!!  false depending on whether aims was linked against ELSI externally or not.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
!!  SOURCE
!!
pure logical function ext_elsi_yes_no() result(y)
  implicit none
  y = .true.
end function ext_elsi_yes_no
!!***
