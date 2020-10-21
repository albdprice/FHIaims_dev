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

AM_CPPFLAGS = -I../sax -I../wxml
AM_FCFLAGS = -I../sax -I../wxml

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
x90_srcs = \
  m_dom_attribute.f90 \
  m_dom_debug.f90 \
  m_dom_document.f90 \
  m_dom_element.f90 \
  m_dom_error.f90 \
  m_dom_namednodemap.f90 \
  m_dom_node.f90 \
  m_dom_nodelist.f90 \
  m_dom_parse.f90 \
  m_dom_types.f90 \
  m_dom_utils.f90 \
  m_dom_strings.f90 \
  xmlf90_dom.f90

# Fortran modules
# Note: never expect they will have the same name as their source
if F90_MOD_UPPERCASE
x90_f03_int_mods = \
  M_DOM_ATTRIBUTE.$(MODEXT) \
  M_DOM_DEBUG.$(MODEXT) \
  M_DOM_DOCUMENT.$(MODEXT) \
  M_DOM_ELEMENT.$(MODEXT) \
  M_DOM_ERROR.$(MODEXT) \
  M_DOM_NAMEDNODEMAP.$(MODEXT) \
  M_DOM_NODE.$(MODEXT) \
  M_DOM_NODELIST.$(MODEXT) \
  M_DOM_PARSE.$(MODEXT) \
  M_DOM_TYPES.$(MODEXT) \
  M_DOM_UTILS.$(MODEXT) \
  M_DOM_STRINGS.$(MODEXT)

x90_f03_ext_mods = \
  XMLF90_DOM.$(MODEXT)
else
x90_f03_int_mods = \
  m_dom_attribute.$(MODEXT) \
  m_dom_debug.$(MODEXT) \
  m_dom_document.$(MODEXT) \
  m_dom_element.$(MODEXT) \
  m_dom_error.$(MODEXT) \
  m_dom_namednodemap.$(MODEXT) \
  m_dom_node.$(MODEXT) \
  m_dom_nodelist.$(MODEXT) \
  m_dom_parse.$(MODEXT) \
  m_dom_types.$(MODEXT) \
  m_dom_utils.$(MODEXT) \
  m_dom_strings.$(MODEXT)

x90_f03_ext_mods = \
  xmlf90_dom.$(MODEXT)
endif

# Libraries to build
noinst_LTLIBRARIES = libxmlf90_dom.la

libxmlf90_dom_la_SOURCES = $(x90_srcs)

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

                    ########################################

# Explicit dependencies within XMLF90
# Note: this is needed because of Fortran
m_dom_attribute.$(LTOBJEXT): \
  m_dom_node.$(LTOBJEXT) \
  m_dom_types.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

m_dom_document.$(LTOBJEXT): \
  m_dom_types.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

m_dom_element.$(LTOBJEXT): \
  m_dom_attribute.$(LTOBJEXT) \
  m_dom_debug.$(LTOBJEXT) \
  m_dom_document.$(LTOBJEXT) \
  m_dom_namednodemap.$(LTOBJEXT) \
  m_dom_node.$(LTOBJEXT) \
  m_dom_nodelist.$(LTOBJEXT) \
  m_dom_types.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

m_dom_namednodemap.$(LTOBJEXT): \
  m_dom_types.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

m_dom_node.$(LTOBJEXT): \
  m_dom_debug.$(LTOBJEXT) \
  m_dom_error.$(LTOBJEXT) \
  m_dom_namednodemap.$(LTOBJEXT) \
  m_dom_nodelist.$(LTOBJEXT) \
  m_dom_types.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

m_dom_nodelist.$(LTOBJEXT): \
  m_dom_types.$(LTOBJEXT)

m_dom_parse.$(LTOBJEXT): \
  m_dom_debug.$(LTOBJEXT) \
  m_dom_document.$(LTOBJEXT) \
  m_dom_element.$(LTOBJEXT) \
  m_dom_node.$(LTOBJEXT) \
  m_dom_types.$(LTOBJEXT)

m_dom_types.$(LTOBJEXT): \
  m_dom_strings.$(LTOBJEXT)

m_dom_utils.$(LTOBJEXT): \
  m_dom_debug.$(LTOBJEXT) \
  m_dom_document.$(LTOBJEXT) \
  m_dom_element.$(LTOBJEXT) \
  m_dom_namednodemap.$(LTOBJEXT) \
  m_dom_node.$(LTOBJEXT) \
  m_dom_types.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

xmlf90_dom.$(LTOBJEXT): \
  m_dom_attribute.$(LTOBJEXT) \
  m_dom_document.$(LTOBJEXT) \
  m_dom_element.$(LTOBJEXT) \
  m_dom_namednodemap.$(LTOBJEXT) \
  m_dom_node.$(LTOBJEXT) \
  m_dom_nodelist.$(LTOBJEXT) \
  m_dom_parse.$(LTOBJEXT) \
  m_dom_types.$(LTOBJEXT) \
  m_dom_utils.$(LTOBJEXT) \
  m_dom_strings.$(LTOBJEXT)

