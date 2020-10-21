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
# External dependencies
#

AM_CPPFLAGS = -I../sax
AM_FCFLAGS = -I../sax

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
x90_srcs = \
  m_xpath_api.f90 \
  xmlf90_xpath.f90

# Other files
EXTRA_DIST = \
  Ba.xml \
  Developer.Guide

# Fortran modules
# Note: never expect they will have the same name as their source
if F90_MOD_UPPERCASE
x90_f03_int_mods = \
  M_XPATH_API.$(MODEXT)

x90_f03_ext_mods = \
  XMLF90_XPATH.$(MODEXT)
else
x90_f03_int_mods = \
  m_xpath_api.$(MODEXT)

x90_f03_ext_mods = \
  xmlf90_xpath.$(MODEXT)
endif

# Libraries to build
noinst_LTLIBRARIES = libxmlf90_xpath.la

libxmlf90_xpath_la_SOURCES = $(x90_srcs)

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(x90_f03_ext_mods) $(x90_f03_int_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(x90_f03_ext_mods) $(x90_f03_int_mods)

# Local cleaning
clean-local:
	rm -f $(x90_f03_int_mods) $(x90_f03_ext_mods) *.obj

                    # ------------------------------------ #

# Basic manual tests
noinst_PROGRAMS = match_tester

match_tester_SOURCES = match_tester.f90
match_tester_LDADD = libxmlf90_xpath.la ../sax/libxmlf90_sax.la

                    # ------------------------------------ #

# Explicit dependencies within XMLF90
# Note: this is needed because of Fortran
xmlf90_xpath.$(LTOBJEXT): \
  m_xpath_api.$(LTOBJEXT)

