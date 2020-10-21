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
# 
#

# Library components, in build order
SUBDIRS = wxml cml sax xpath dom .

# Libraries to install
lib_LTLIBRARIES = libxmlf90.la

nodist_libxmlf90_la_SOURCES = xmlf90_info.f90
libxmlf90_la_LDFLAGS = -version-info 0:0:0
libxmlf90_la_LIBADD = \
  cml/libxmlf90_cml.la \
  dom/libxmlf90_dom.la \
  sax/libxmlf90_sax.la \
  wxml/libxmlf90_wxml.la \
  xpath/libxmlf90_xpath.la

# Fortran modules
# Note: never expect they will have the same name as their source
if F90_MOD_UPPERCASE
x90_f03_mods = \
  XMLF90_INFO.$(MODEXT)
else
x90_f03_mods = \
  xmlf90_info.$(MODEXT)
endif

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(x90_f03_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(x90_f03_mods)

# Local cleaning
CLEANFILES = $(x90_f03_mods)
