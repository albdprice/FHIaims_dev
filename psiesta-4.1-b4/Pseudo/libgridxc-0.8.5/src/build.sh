#!/bin/sh
#
# This script harnesses the standard libGridXC makefile to build an
# alternative hierarchy of the form
#
#                PREFIX/
#                      gridxc.mk
#                      libxc.mk
#                      serial/
#                            gridxc.mk
#                            libxc.mk
#                            lib/
#                            include/
#                      mpi/
#                            gridxc.mk
#                            libxc.mk
#                            lib/
#                            include/
#
# The 'mpi' section is optional, depending on the command issued.
#
# Usage:
#
# WITH_LIBXC=X  WITH_MPI=Y PREFIX=Z sh build.sh
#
# If Y is empty, only the 'serial' section is built, otherwise both
# 'serial' and 'mpi' sections are created. If X is non-empty, libXC
# support is compiled in. Note that there are no extra sections
# devoted to libXC/MPI combinations yet. The user should use separate
# top-levels for this.
#
# If PREFIX is empty, the hierarchy is rooted in the current directory.
# 
# The top-level gridxc.mk has special logic to dispatch the
# appropriate variables depending on the setting of WITH_MPI in client
# makefiles.
#
# The individual sections 'serial' and 'mpi' are still bona-fide
# libGridXC installations in the old style.
#
#
# Set to "." if empty or unset
inPREFIX=${PREFIX:-.}
#
# Build first without MPI
#
echo "==> make WITH_LIBXC=${WITH_LIBXC} WITH_MPI= PREFIX=${inPREFIX}/serial"
sleep 1
make clean
make WITH_LIBXC=${WITH_LIBXC} WITH_MPI= PREFIX=${inPREFIX}/serial
#
# Install the top-level gridxc.mk
#
if [ "${WITH_LIBXC}" != "" ]
then
    sed 's/GRIDXC_USES_LIBXC=0/GRIDXC_USES_LIBXC=1/g' top.gridxc.mk.in \
	                    > ${inPREFIX}/gridxc.mk
else
    cp -p top.gridxc.mk.in ${inPREFIX}/gridxc.mk
fi
#
# Build with MPI if requested
#
if [ "${WITH_MPI}" != "" ]
then
   echo "==> make WITH_LIBXC=${WITH_LIBXC} WITH_MPI=1 PREFIX=${inPREFIX}/mpi"
   sleep 1
   make clean
   make WITH_LIBXC=${WITH_LIBXC} WITH_MPI=1 PREFIX=${inPREFIX}/mpi
#
# Install the top-level gridxc.mk
#
   if [ "${WITH_LIBXC}" != "" ]
   then
       sed 's/GRIDXC_USES_LIBXC=0/GRIDXC_USES_LIBXC=1/g' top.gridxc.mk.in | \
	   sed 's/GRIDXC_USES_MPI=0/GRIDXC_USES_MPI=1/g'  > ${inPREFIX}/gridxc.mk
   else
       sed 's/GRIDXC_USES_MPI=0/GRIDXC_USES_MPI=1/g' top.gridxc.mk.in \
	                         > ${inPREFIX}/gridxc.mk
   fi
fi
#
#   Install copy of libxc.mk file at the top
#
cp -p libxc.mk ${inPREFIX}/libxc.mk


