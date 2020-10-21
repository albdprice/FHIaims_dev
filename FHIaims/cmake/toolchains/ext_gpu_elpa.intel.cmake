set(TARGET_NAME aims.x CACHE STRING "")
set(CMAKE_INSTALL_PREFIX "$ENV{AIMS_HOME}" CACHE STRING "")

###############
# Basic Flags #
###############
set(CMAKE_Fortran_COMPILER mpiifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(Fortran_MIN_FLAGS "-O0 -fp-model precise" CACHE STRING "")
set(INC_PATHS "$ENV{ELPAROOT}/include/elpa-2019.05.001/modules" CACHE STRING "")
set(LIB_PATHS "$ENV{ELPAROOT}/lib $ENV{MKLROOT}/lib/intel64 $ENV{CUDA_HOME}/lib64" CACHE STRING "")
set(LIBS "elpa mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_sequential mkl_core cublas cudart" CACHE STRING "")

###############
# C,C++ Flags #
###############
set(USE_CXX_FILES ON CACHE BOOL "")
set(CMAKE_C_COMPILER icc CACHE STRING "")
set(CMAKE_C_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")
set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
set(CMAKE_CXX_FLAGS "-O3 -ip -fp-model precise" CACHE STRING "")

##############
# ELPA Flags #
##############
set(USE_EXTERNAL_ELPA ON CACHE BOOL "")
set(USE_GPU_ELPA ON CACHE BOOL "")

##########################
# GPU Acceleration Flags #
##########################
set(USE_CUDA ON CACHE BOOL "")
set(CMAKE_CUDA_COMPILER nvcc CACHE STRING "")
set(CMAKE_CUDA_FLAGS "-O3 -m64 -DAdd_ -arch=sm_60 -lcublas" CACHE STRING "")
set(CMAKE_CUDA_HOST_COMPILER g++ CACHE STRING "")
