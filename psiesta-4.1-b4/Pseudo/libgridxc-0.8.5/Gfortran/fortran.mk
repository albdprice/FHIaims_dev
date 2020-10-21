#
# Example Fortran macros: gfortran
#
# Make sure this variable is set if you intend to use libxc
# You can set it here, or elsewhere (e.g., in an environment module)
### LIBXC_ROOT=/path/to/libxc
#
# These two instances are needed, instead of just FC
FC_SERIAL=gfortran
FC_PARALLEL=mpif90
#

FFLAGS= -O2  -fimplicit-none
FFLAGS_DEBUG= -g -O0 -fbacktrace -fcheck=all -fimplicit-none
FFLAGS_CHECKS= -g  -O0 -g -Wall -Wextra -Warray-temporaries \
              -Wconversion -fimplicit-none -fbacktrace \
              -ffree-line-length-0 -fcheck=all \
              -ffpe-trap=zero,overflow,underflow -finit-real=nan
LDFLAGS=
#
AR=ar
RANLIB=ranlib
#
DEFS_PREFIX=
INC_PREFIX= -I
MOD_PREFIX= -I
MOD_EXT=.mod
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#










