!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Print the build configuration of the current instance of FHI-aims.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
subroutine cmake_info()

  use localorb_io, only: localorb_multi

  implicit none

  character(:), allocatable :: string, word
  logical :: link_first_line

  link_first_line = .true.

  call localorb_multi( &
       & 'Build configuration of the current instance of FHI-aims', &
       & '-------------------------------------------------------', &
       & 'FHI-aims version      : 190903', &
       & 'Commit number         : 82356777a', &
       & 'CMake host system     : Linux-4.15.0-1059-oem', &
       & 'CMake version         : 3.10.2', &
       & 'Fortran compiler      : /opt/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpif90 &
       &(Intel) version &
       &19.0.5.20190815', &
       & 'Fortran compiler flags: -O3 -ip -fp-model precise', &
       format='(2x, a)')
  if ('/opt/intel/bin/icc' /= '') &
       & call localorb_multi( &
       & 'C compiler            : /opt/intel/bin/icc &
       &(Intel) version 19.0.5.20190815', &
       & 'C compiler flags      : -O3 -ip -fp-model precise', &
       & format='(2x, a)')
  if ('' /= '') &
       & call localorb_multi( &
       & 'C++ compiler          :  &
       &() version ', &
       & 'C++ compiler flags    : ', &
       & format='(2x, a)')
  if ('' /= '') &
       & call localorb_multi( &
       & 'CUDA compiler         : ', &
       & 'CUDA compiler flags   : ', &
       & format='(2x, a)')
  call localorb_multi('Architecture          : ', &
       & format='(2x, a)')
  if ('ON' == 'ON') &
       & call localorb_multi('Using MPI', format='(2x, a)')
  if ('ON' == 'ON') &
       & call localorb_multi('Using Scalapack', format='(2x, a)')
  if ('ON' == 'ON') &
       & call localorb_multi('Using C files', format='(2x, a)')
  if ('ON' == 'ON') &
       & call localorb_multi('Using LibXC', format='(2x, a)')
  if ('OFF' == 'ON') &
       & call localorb_multi('Using CUDA', format='(2x, a)')
  if ('OFF' == 'ON') &
       & call localorb_multi('Using CFFI', format='(2x, a)')
  if ('ON' == 'ON') &
       & call localorb_multi('Using SPGlib', format='(2x, a)')
  if ('ON' == 'ON') &
       & call localorb_multi('Using i-PI', format='(2x, a)')
  if ('OFF' == 'ON') &
       & call localorb_multi('Using HDF5', format='(2x, a)')
  if ('' == 'ON') &
       & call localorb_multi('Using external ELSI', format='(2x, a)')
  if ('' == 'ON') &
       & call localorb_multi('Using external ELPA', format='(2x, a)')
  if ('OFF' == 'ON') &
       & call localorb_multi('Using GPU ELPA as default', format='(2x, a)')
  ! Print a list of linked libraries
  string = '/opt/intel/mkl/lib/intel64/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64/libmkl_sequential.so;/opt/intel/mkl/lib/intel64/libmkl_core.so;/opt/intel/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.so;/opt/intel/mkl/lib/intel64/libmkl_scalapack_lp64.so;'
  do while (index(string,';') > 0)
     word = string(:index(string,';')-1)
     if (link_first_line) then
        call localorb_multi('Linking against: '//word, format='(2x, a)')
        link_first_line = .false.
     else
        call localorb_multi(word, format='(19x, a)')
     end if
     string = string(index(string,';')+1:)
  end do
  call localorb_multi('')
end subroutine cmake_info
