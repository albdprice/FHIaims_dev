#
# Makefile for libGridXC 
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
#
# -------- This statement is important
#          for the operation of sub-makes
#   It works with GNUmake
#
export
#----------------------------------------
# Clean first in all cases, to avoid
# inconsistent LibXC/mpi states
#
default:  clean module
#
VPATH=.
MAIN_OBJDIR=.
#
#
PREFIX ?= $(MAIN_OBJDIR)
GRIDXC_PREFIX=$(PREFIX)
#
#
.PHONY: module setup_mk_file
#
# We allow for an external specification of the location 
# of the arch.make file. This must be an absolute path
#
#
ARCH_MAKE_DEFAULT=$(MAIN_OBJDIR)/fortran.mk
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)

# Ensure that we can ALWAYS run make dep in the MAIN directory
ifeq ($(VPATH),.)
 -include $(ARCH_MAKE)
else
 include $(ARCH_MAKE)
endif

# Default library commands
AR ?= ar
RANLIB ?= ranlib

ifdef WITH_LIBXC

 ifndef LIBXC_ROOT
   $(error you need to define LIBXC_ROOT in your arch.make)
 endif

 FPPFLAGS += -DLIBXC
 include libxc.mk

endif

dep:
	sfmakedepend --depend=obj  --modext=o *.f *.f90 *.F *.F90

GRID_SRCS= alloc.F90 bessph.f cellsubs.f \
             chkgmx.f  interpolation.f90  \
             minvec.f m_io.f moreParallelSubs.F90 \
             m_fft_gpfa.F  m_walltime.f90 \
	     precision.F radfft.f sorting.f sys.F \
             array.F90
GRID_OBJS:= $(GRID_SRCS:.f=.o)
GRID_OBJS:= $(GRID_OBJS:.F=.o)
GRID_OBJS:= $(GRID_OBJS:.f90=.o)
GRID_OBJS:= $(GRID_OBJS:.F90=.o)

LOCAL_SRCS= atomxc.F90 cellxc.F90 debugxc.f90 ggaxc.F ldaxc.F \
            gridxc.F90 vdwxc.F90 xcmod.f90  fft3d.F90 fftr.F90 \
            mesh1d.f90 mesh3d.F90 am05.f90 vv_vdwxc.F90 \
            gridxc_config.F90
LOCAL_OBJS:= $(LOCAL_SRCS:.f=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.F=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.f90=.o)
LOCAL_OBJS:= $(LOCAL_OBJS:.F90=.o)

ALL_OBJS=$(GRID_OBJS) $(LOCAL_OBJS)


INCFLAGS:= $(LIBXC_INCFLAGS) $(INCFLAGS)

ifdef WITH_MPI
  FC=$(FC_PARALLEL)
  MPI_INTERFACE=libmpi_f90.a
  FPPFLAGS:= -DMPI $(FPPFLAGS)
  $(ALL_OBJS):  INCFLAGS += -I./MPI_instr
else
  FC=$(FC_SERIAL)
  MPI_INTERFACE=
endif

libmpi_f90.a: 
	@(cd MPI_instr ; $(MAKE) \
                        "VPATH=$(VPATH)/MPI_instr" )

libGridXC.a: $(MPI_INTERFACE) $(ALL_OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libGridXC.a $(ALL_OBJS) 
	@if [ ! -z "$(MPI_INTERFACE)" ] ; then \
            $(AR) $(ARFLAGS_EXTRA) cru libGridXC.a MPI_instr/*.o ; fi
	-$(RANLIB) libGridXC.a

MODULES_TO_INSTALL=gridxc.mod gridxc_config.mod m_atomxc.mod \
                   m_cellxc.mod xcmod.mod m_ldaxc.mod m_ggaxc.mod \
                   m_fft_gpfa.mod debugxc.mod mesh3d.mod

ifdef WITH_LIBXC
 ifdef WITH_MPI
setup_mk_file:
	sed 's/GRIDXC_USES_LIBXC=0/GRIDXC_USES_LIBXC=1/g' gridxc.mk.in | \
	sed 's/GRIDXC_USES_MPI=0/GRIDXC_USES_MPI=1/g'  > gridxc.mk
 else	
setup_mk_file:
	sed 's/GRIDXC_USES_LIBXC=0/GRIDXC_USES_LIBXC=1/g' gridxc.mk.in > gridxc.mk
 endif
else
 ifdef WITH_MPI
setup_mk_file:
	sed 's/GRIDXC_USES_MPI=0/GRIDXC_USES_MPI=1/g' gridxc.mk.in > gridxc.mk
 else	
setup_mk_file:
	cp gridxc.mk.in gridxc.mk
 endif
endif

module: setup_mk_file libGridXC.a
	@echo "Prefix: $(PREFIX)"
	mkdir -p lib 
	mkdir -p include
	cp -p $(MODULES_TO_INSTALL) include/
	cp -p libGridXC.a lib/
	@if [ "$(PREFIX)" != "$(MAIN_OBJDIR)" ]; then \
          echo "mkdir -p $(GRIDXC_PREFIX)/lib" ; \
          mkdir -p $(GRIDXC_PREFIX)/lib ; \
          echo "mkdir -p $(GRIDXC_PREFIX)/include" ; \
          mkdir -p $(GRIDXC_PREFIX)/include ; \
          echo "cp -p $(MODULES_TO_INSTALL) $(GRIDXC_PREFIX)/include/" ; \
          cp -p $(MODULES_TO_INSTALL) $(GRIDXC_PREFIX)/include/ ; \
          echo "cp -p libGridXC.a $(GRIDXC_PREFIX)/lib" ; \
          cp -p libGridXC.a $(GRIDXC_PREFIX)/lib ; \
          echo "cp -p gridxc.mk libxc.mk $(GRIDXC_PREFIX)" ; \
          cp -p gridxc.mk libxc.mk $(GRIDXC_PREFIX) ; fi

clean:
	rm -f *.o *.*d *.a gridxc.mk
	rm -rf lib include
	(cd MPI_instr ; $(MAKE) \
                        "VPATH=$(VPATH)/MPI_instr" clean)


#
# DO NOT DELETE THIS LINE - used by make depend
am05.o: precision.o sys.o
array.o: precision.o sys.o
atomxc.o: alloc.o debugxc.o ggaxc.o gridxc_config.o ldaxc.o mesh1d.o
atomxc.o: precision.o radfft.o sys.o vdwxc.o xcmod.o
bessph.o: precision.o sys.o
cellsubs.o: precision.o
cellxc.o: alloc.o cellsubs.o chkgmx.o debugxc.o fftr.o ggaxc.o gridxc_config.o
cellxc.o: ldaxc.o mesh3d.o moreParallelSubs.o precision.o sys.o vdwxc.o xcmod.o
chkgmx.o: cellsubs.o minvec.o precision.o sys.o
debugxc.o: m_io.o moreParallelSubs.o
fft3d.o: alloc.o debugxc.o m_fft_gpfa.o mesh3d.o precision.o sys.o
fftr.o: alloc.o fft3d.o mesh3d.o precision.o sys.o
ggaxc.o: am05.o ldaxc.o precision.o sys.o
gridxc.o: atomxc.o cellxc.o debugxc.o ggaxc.o gridxc_config.o ldaxc.o
gridxc.o: m_fft_gpfa.o mesh3d.o precision.o xcmod.o
gridxc_config.o: sys.o
ldaxc.o: precision.o sys.o
m_io.o: sys.o
mesh1d.o: interpolation.o precision.o
mesh3d.o: alloc.o array.o debugxc.o gridxc_config.o precision.o sorting.o sys.o
minvec.o: cellsubs.o precision.o sorting.o sys.o
moreParallelSubs.o: alloc.o array.o gridxc_config.o m_io.o precision.o sys.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
vdwxc.o: alloc.o debugxc.o ggaxc.o interpolation.o ldaxc.o mesh1d.o precision.o
vdwxc.o: radfft.o sys.o vv_vdwxc.o
vv_vdwxc.o: alloc.o debugxc.o interpolation.o mesh1d.o precision.o radfft.o
vv_vdwxc.o: sys.o
xcmod.o: precision.o sys.o vdwxc.o
m_array.o: array.o
m_atomxc.o: atomxc.o
m_bessph.o: bessph.o
m_cellxc.o: cellxc.o
m_chkgmx.o: chkgmx.o
m_fft3d.o: fft3d.o
m_ggaxc.o: ggaxc.o
m_ldaxc.o: ldaxc.o
m_minvec.o: minvec.o
m_radfft.o: radfft.o
m_vdwxc.o: vdwxc.o
m_vv_vdwxc.o: vv_vdwxc.o
