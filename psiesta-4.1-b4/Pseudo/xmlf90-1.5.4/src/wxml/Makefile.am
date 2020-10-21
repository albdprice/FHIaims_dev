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
# Main source files
#

# Common source files
x90_srcs = \
  m_wxml_array_str.f90 \
  m_wxml_buffer.f90 \
  m_wxml_core.f90 \
  m_wxml_dictionary.f90 \
  m_wxml_elstack.f90 \
  m_wxml_error.f90 \
  m_wxml_escape.f90 \
  m_wxml_overloads.f90 \
  m_wxml_text.F90 \
  xmlf90_wxml.f90

# Other files
EXTRA_DIST = \
  README

# Fortran modules
# Note: never expect they will have the same name as their source
if F90_MOD_UPPERCASE
x90_f03_int_mods = \
  M_WXML_ARRAY_STR.$(MODEXT) \
  M_WXML_BUFFER.$(MODEXT) \
  M_WXML_CORE.$(MODEXT) \
  M_WXML_DICTIONARY.$(MODEXT) \
  M_WXML_ELSTACK.$(MODEXT) \
  M_WXML_ERROR.$(MODEXT) \
  M_WXML_ESCAPE.$(MODEXT) \
  M_WXML_OVERLOADS.$(MODEXT) \
  M_WXML_TEXT.$(MODEXT)

x90_f03_ext_mods = \
  XMLF90_WXML.$(MODEXT)
else
x90_f03_int_mods = \
  m_wxml_array_str.$(MODEXT) \
  m_wxml_buffer.$(MODEXT) \
  m_wxml_core.$(MODEXT) \
  m_wxml_dictionary.$(MODEXT) \
  m_wxml_elstack.$(MODEXT) \
  m_wxml_error.$(MODEXT) \
  m_wxml_escape.$(MODEXT) \
  m_wxml_overloads.$(MODEXT) \
  m_wxml_text.$(MODEXT)

x90_f03_ext_mods = \
  xmlf90_wxml.$(MODEXT)
endif

# Libraries to build
noinst_LTLIBRARIES = libxmlf90_wxml.la

libxmlf90_wxml_la_SOURCES = $(x90_srcs)

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
m_wxml_buffer.$(LTOBJEXT): \
  m_wxml_error.$(LTOBJEXT)

m_wxml_core.$(LTOBJEXT): \
  m_wxml_buffer.$(LTOBJEXT) \
  m_wxml_dictionary.$(LTOBJEXT) \
  m_wxml_elstack.$(LTOBJEXT)

m_wxml_dictionary.$(LTOBJEXT): \
  m_wxml_array_str.$(LTOBJEXT) \
  m_wxml_escape.$(LTOBJEXT)

m_wxml_elstack.$(LTOBJEXT): \
  m_wxml_error.$(LTOBJEXT)

m_wxml_overloads.$(LTOBJEXT): \
  m_wxml_core.$(LTOBJEXT) \
  m_wxml_text.$(LTOBJEXT)

xmlf90_wxml.$(LTOBJEXT): \
  m_wxml_core.$(LTOBJEXT) \
  m_wxml_overloads.$(LTOBJEXT) \
  m_wxml_text.$(LTOBJEXT)

