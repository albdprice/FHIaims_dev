! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Return details about ELSI's versioning.
!!
subroutine elsi_version_info(version,datestamp,commit,hostname,datetime)

   implicit none

   character(len=8), intent(out) :: version
   character(len=8), intent(out) :: datestamp
   character(len=8), intent(out) :: commit
   character(len=40), intent(out) :: hostname
   character(len=20), intent(out) :: datetime

   version = "2.3.1"
   datestamp = "20190714"
   commit = "82356777a"
   hostname = "Linux-4.15.0-1059-oem"
   datetime = "2019-11-05T20:15:49Z"

end subroutine

!>
!! Return details about ELSI's installation.
!!
subroutine elsi_solver_info(have_aeo,have_pexsi,have_sips)

   implicit none

   logical, intent(out) :: have_aeo
   logical, intent(out) :: have_pexsi
   logical, intent(out) :: have_sips

   character(len=*), parameter :: tmp1 = "OFF"
   character(len=*), parameter :: tmp2 = "OFF"
   character(len=*), parameter :: tmp3 = "OFF"

   if(tmp1 == "ON") then
      have_aeo = .true.
   else
      have_aeo = .false.
   end if

   if(tmp2 == "ON") then
      have_pexsi = .true.
   else
      have_pexsi = .false.
   end if

   if(tmp3 == "ON") then
      have_sips = .true.
   else
      have_sips = .false.
   end if

end subroutine
