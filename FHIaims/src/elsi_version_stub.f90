! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.
!

!>
!! This subroutine returns details about ELSI's versioning.
!!
subroutine elsi_version_info(version,datestamp,commit,hostname,datetime)

   implicit none

   character(len=8), intent(out) :: version
   character(len=8), intent(out) :: datestamp
   character(len=8), intent(out) :: commit
   character(len=40), intent(out) :: hostname
   character(len=20), intent(out) :: datetime

   version = "FHI-aims"
   datestamp = "N/A"
   commit = "N/A"
   hostname = "N/A"
   datetime = "N/A"

end subroutine

!>
!! This subroutine returns details about ELSI's solver libraries.
!!
subroutine elsi_solver_info(have_aeo,have_pexsi,have_sips)

   implicit none

   logical, intent(out) :: have_aeo
   logical, intent(out) :: have_pexsi
   logical, intent(out) :: have_sips

   have_aeo = .false.
   have_pexsi = .false.
   have_sips = .false.

end subroutine
