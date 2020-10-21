#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#
GRIDXC_ROOT=..
#
# Get building rules for this compiler
#
include ../fortran.mk
#
# Find appropriate modules and libraries
#
include ../gridxc.mk
#
ifeq ($(GRIDXC_USES_LIBXC),1)
default: test1 test2 test3 test4 test-libxc test-libxc2 test-match
else
default: test1 test2 test3 test4
endif
#
ifeq ($(GRIDXC_USES_MPI),1)
  $(info  -- Will use by default FC=mpif90. You can override)
  FC = mpif90
  FPPFLAGS := -DMPI $(FPPFLAGS)
endif

#
# The GRIDXC_ variables are set by the gridxc.mk file
#
INCFLAGS := $(GRIDXC_INCFLAGS) $(INCFLAGS)
LIBS:= $(GRIDXC_LIBS) $(LIBS)

#
test1:  test1.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test1 test1.o handlers.o $(LIBS)
#
test2:  test2.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test2 test2.o handlers.o $(LIBS)
#
test3:  test3.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test3 test3.o handlers.o $(LIBS)
#
test4:  test4.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test4 test4.o handlers.o $(LIBS)
#
test-match:  test-match.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test-match test-match.o handlers.o $(LIBS)
#
test-libxc:  test-libxc.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test-libxc test-libxc.o handlers.o $(LIBS)
test-libxc2:  test-libxc2.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test-libxc2 test-libxc2.o handlers.o $(LIBS)
#
clean:
	rm -f *.o *.*d *.a test1 test2 test3 test4 test-match
	rm -f test-libxc test-libxc2
	rm -f atomVxc.out cellVxc.out  debugXC*out
	rm -f vdw_kernel.table*  *.Vxc


