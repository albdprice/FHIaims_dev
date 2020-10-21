# GRIDXC helper mk
#
# In order to hide the structure of this hierarchy, this file
# will export the needed symbols for a client.
#
# All you need is to define GRIDXC_ROOT in your makefile
# to point to your installation of the gridxc library.
#
#----------------------------------------------------------
#
GRIDXC_USES_LIBXC=0
GRIDXC_USES_MPI=0
#
ifeq ($(GRIDXC_USES_MPI),1)
  $(info GRIDXC was compiled with MPI)
endif
#
ifeq ($(GRIDXC_USES_LIBXC),1)
 $(info GRIDXC was compiled with libxc support)
 #
 ifndef LIBXC_ROOT
   $(error you need to define LIBXC_ROOT in your arch.make)
 endif
 #
 include $(GRIDXC_ROOT)/libxc.mk
 #
else
 # 
 LIBXC_INCFLAGS=
 LIBXC_LIBS=
endif
#
GRIDXC_INCFLAGS= -I $(GRIDXC_ROOT)/include $(LIBXC_INCFLAGS)
GRIDXC_LIBS=$(GRIDXC_ROOT)/lib/libGridXC.a $(LIBXC_LIBS)
