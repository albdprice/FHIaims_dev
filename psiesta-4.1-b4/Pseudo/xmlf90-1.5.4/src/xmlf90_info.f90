module xmlf90_info

  implicit none

  private

  public :: xmlf90_info_version

  integer, parameter :: XMLF90_VERSION_MAJOR = 1
  integer, parameter :: XMLF90_VERSION_MINOR = 5
  integer, parameter :: XMLF90_VERSION_MICRO = 4

contains

  !> @brief Provides the version number triplet of XMLF90
  !!
  !! @param[out] major: major version number
  !! @param[out] minor: minor version number
  !! @param[out] micro: micro version number
  subroutine xmlf90_info_version(major, minor, micro)

    implicit none

    integer, intent(out) :: major, minor, micro

    major = XMLF90_VERSION_MAJOR
    minor = XMLF90_VERSION_MINOR
    micro = XMLF90_VERSION_MICRO

  end subroutine xmlf90_info_version

end module xmlf90_info
