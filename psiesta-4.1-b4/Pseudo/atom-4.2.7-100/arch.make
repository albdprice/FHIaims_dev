#
# Sample arch.make for stand-alone ATOM program
#
ARCH=unknown
#
# --- Required libraries: xmlf90 and libGridXC
#
# It is assumed that the user has downloaded the library tarballs from Launchpad
# and compiled and installed them.
#
# Somehow (here, or in the environment) define the XMLF90_ROOT and
# GRIDXC_ROOT symbols
#
XMLF90_ROOT=/home/lleblanc/src/psiesta-4.1-b4/Pseudo/xmlf90-1.5.4/
GRIDXC_ROOT=/home/lleblanc/src/psiesta-4.1-b4/Pseudo/libgridxc-0.8.5/Gfortran
#
#  NOTE: The building mechanism for Siesta-related libraries is still
#  being refined. The paths of the .mk files under $(XMLF90_ROOT) and $(GRIDXC_ROOT)
#  might change with different versions of the libraries. Check their
#  installation trees and use the appropriate path to the .mk file.
#
include $(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk
include $(GRIDXC_ROOT)/gridxc.mk
#
# -----Compiler-dependent settings -------------------------------------
#
FC=gfortran
#
FFLAGS= -O2
FFLAGS_DEBUG= -O0 -g -fbacktrace
LDFLAGS=
RANLIB=echo
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
