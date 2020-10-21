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

AM_CPPFLAGS = -I../wxml
AM_FCFLAGS = -I../wxml

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
x90_srcs = \
  m_wcml.f90 \
  m_wstml.f90 \
  xmlf90_cml.f90

# Fortran modules
# Note: never expect they will have the same name as their source
if F90_MOD_UPPERCASE
x90_f03_int_mods = \
  M_WCML.$(MODEXT) \
  M_WSTML.$(MODEXT)

x90_f03_ext_mods = \
  XMLF90_CML.$(MODEXT)
else
x90_f03_int_mods = \
  m_wcml.$(MODEXT) \
  m_wstml.$(MODEXT)

x90_f03_ext_mods = \
  xmlf90_cml.$(MODEXT)
endif

# Libraries to build
noinst_LTLIBRARIES = libxmlf90_cml.la

libxmlf90_cml_la_SOURCES = $(x90_srcs)

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

# Explicit dependencies within XMLF90
# Note: this is needed because of Fortran
m_wcml.$(LTOBJEXT): \
  m_wstml.$(LTOBJEXT)

xmlf90_cml.$(LTOBJEXT): \
  m_wcml.$(LTOBJEXT) \
  m_wstml.$(LTOBJEXT)

