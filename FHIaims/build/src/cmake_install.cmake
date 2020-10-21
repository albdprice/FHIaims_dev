# Install script for directory: /home/albd/research/FHIaims/running_copy/FHIaims/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/albd/research/FHIaims/running_copy/FHIaims/build/aims.190903.scalapack.mpi.x")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x"
         OLD_RPATH "/opt/intel/mkl/lib/intel64:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aims.190903.scalapack.mpi.x")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/aims/aimsConfig.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/aims/aimsConfig.cmake"
         "/home/albd/research/FHIaims/running_copy/FHIaims/build/src/CMakeFiles/Export/lib/cmake/aims/aimsConfig.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/aims/aimsConfig-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/aims/aimsConfig.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/aims" TYPE FILE FILES "/home/albd/research/FHIaims/running_copy/FHIaims/build/src/CMakeFiles/Export/lib/cmake/aims/aimsConfig.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/aims" TYPE FILE FILES "/home/albd/research/FHIaims/running_copy/FHIaims/build/src/CMakeFiles/Export/lib/cmake/aims/aimsConfig-noconfig.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/aims" TYPE FILE FILES "/home/albd/research/FHIaims/running_copy/FHIaims/build/src/aimsConfigVersion.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/external/elsi/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/subdirectory_aims1/build/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/subdirectory_aims2/build/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/subdirectory_aims3/build/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/subdirectory_aims4/build/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/subdirectory_aims5/build/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/basis_sets/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/ipc/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/elpa/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/optical_response/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/MagneticResponse/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/soc/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/SCGW/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/Periodic_PT2/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/LRC_PT2/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/mbd-std/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/mbd-dev/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DMFT_embed/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT_reduce_memory/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT_polarizability/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT_phonon/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT_phonon_reduce_memory/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT_dielectric/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/DFPT_phonon_gamma/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/sBGE2/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/RRS-PBC/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/xc_dfauto/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/FCIQMC/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/bse/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/CCcalc/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/qpe_osPT2/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/friction/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/external/cmake_install.cmake")
  include("/home/albd/research/FHIaims/running_copy/FHIaims/build/src/relativity/cmake_install.cmake")

endif()

