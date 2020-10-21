set(TARGET_NAME aims.x CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

###############
# Basic Flags #
###############
set(CMAKE_Fortran_COMPILER $ENV{MPI_ROOT}/bin/mpifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0 -fp-model precise" CACHE STRING "")
set(USE_SCALAPACK OFF CACHE BOOL "")
set(LIB_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "")
set(LIBS "mkl_intel_lp64 mkl_sequential mkl_core" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(CMAKE_C_COMPILER icc CACHE STRING "")
set(CMAKE_C_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")