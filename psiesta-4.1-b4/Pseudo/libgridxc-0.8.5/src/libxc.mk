# libxc helper mk
#
# In order to hide the structure of the libxc installation hierarchy, this file
# will export the needed symbols for a client.
#
# All you need is to define LIBXC_ROOT in your makefile
# or elsewhere (e.g., in an environment module)
# to point to your installation of the libxc library.
#
# The symbols below are valid for libxc 2.2.0-->4.2.3, at least.
# Update this file accordingly if the details change.
# In particular, check the library names.
#
LIBXC_INCFLAGS= -I $(LIBXC_ROOT)/include
LIBXC_LIBS=-L$(LIBXC_ROOT)/lib -l xcf90 -l xc
#
# The line below might be needed in some systems
#
#LIBXC_LIBS=-L$(LIBXC_ROOT)/lib -Wl,-rpath=$(LIBXC_ROOT)/lib -l xcf90 -l xc
