# -*- Automake -*-
#
# Makefile for the XMLF90 package
#
# Copyright (C) 2017 Yann Pouillon
#
# This file is part of the XMLF90 software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Top makefile
#

# Autotools parameters
ACLOCAL_AMFLAGS = -I config/m4

# Subdirectories containing makefiles
SUBDIRS = src doc

# Support for pkg-config (see http://pkg-config.freedesktop.org/wiki/)
pkgconfigdir = $(prefix)/lib/pkgconfig
pkgconfig_DATA = config/data/xmlf90.pc

# Support for SIESTA's .mk files
siestadir = $(datarootdir)/org.siesta-project
siesta_DATA = config/data/xmlf90.mk

# Files to distribute in the source package
EXTRA_DIST = \
  COPYING \
  KNOWN_ISSUES \
  LICENSE \
  README.md \
  config/data/xmlf90.mk.in \
  config/data/xmlf90.pc.in \
  config/data/xmlf90-config.yml.in

# Custom targets
distclean-local:
	rm -f $(distdir).tar.bz2 $(distdir).tar.gz $(distdir).tar.xz
