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
  m_sax_buffer.f90 \
  m_sax_charset.f90 \
  m_sax_converters.f90 \
  m_sax_debug.f90 \
  m_sax_dictionary.f90 \
  m_sax_elstack.f90 \
  m_sax_entities.f90 \
  m_sax_fsm.f90 \
  m_sax_reader.f90 \
  m_xml_error.f90 \
  m_xml_parser.f90 \
  xmlf90_sax.f90

# Fortran modules
# Note: never expect they will have the same name as their source
if F90_MOD_UPPERCASE
x90_f03_int_mods = \
  M_SAX_BUFFER.$(MODEXT) \
  M_SAX_CHARSET.$(MODEXT) \
  M_SAX_CONVERTERS.$(MODEXT) \
  M_SAX_DEBUG.$(MODEXT) \
  M_SAX_DICTIONARY.$(MODEXT) \
  M_SAX_ELSTACK.$(MODEXT) \
  M_SAX_ENTITIES.$(MODEXT) \
  M_SAX_FSM.$(MODEXT) \
  M_SAX_READER.$(MODEXT) \
  M_XML_ERROR.$(MODEXT) \
  M_XML_PARSER.$(MODEXT)

x90_f03_ext_mods = \
  XMLF90_SAX.$(MODEXT)
else
x90_f03_int_mods = \
  m_sax_buffer.$(MODEXT) \
  m_sax_charset.$(MODEXT) \
  m_sax_converters.$(MODEXT) \
  m_sax_debug.$(MODEXT) \
  m_sax_dictionary.$(MODEXT) \
  m_sax_elstack.$(MODEXT) \
  m_sax_entities.$(MODEXT) \
  m_sax_fsm.$(MODEXT) \
  m_sax_reader.$(MODEXT) \
  m_xml_error.$(MODEXT) \
  m_xml_parser.$(MODEXT)

x90_f03_ext_mods = \
  xmlf90_sax.$(MODEXT)
endif

# Libraries to build
noinst_LTLIBRARIES = libxmlf90_sax.la

libxmlf90_sax_la_SOURCES = $(x90_srcs)

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
m_sax_converters.$(LTOBJEXT): \
  m_sax_debug.$(LTOBJEXT)

m_sax_dictionary.$(LTOBJEXT): \
  m_sax_buffer.$(LTOBJEXT)

m_sax_elstack.$(LTOBJEXT): \
  m_sax_buffer.$(LTOBJEXT)

m_sax_entities.$(LTOBJEXT): \
  m_sax_buffer.$(LTOBJEXT) \
  m_xml_error.$(LTOBJEXT)

m_sax_fsm.$(LTOBJEXT): \
  m_sax_buffer.$(LTOBJEXT) \
  m_sax_charset.$(LTOBJEXT) \
  m_sax_dictionary.$(LTOBJEXT) \
  m_sax_elstack.$(LTOBJEXT) \
  m_sax_entities.$(LTOBJEXT)

m_xml_error.$(LTOBJEXT): \
  m_sax_elstack.$(LTOBJEXT)

m_xml_parser.$(LTOBJEXT): \
  m_sax_buffer.$(LTOBJEXT) \
  m_sax_debug.$(LTOBJEXT) \
  m_sax_dictionary.$(LTOBJEXT) \
  m_sax_elstack.$(LTOBJEXT) \
  m_sax_fsm.$(LTOBJEXT) \
  m_sax_reader.$(LTOBJEXT) \
  m_xml_error.$(LTOBJEXT)

xmlf90_sax.$(LTOBJEXT): \
  m_sax_converters.$(LTOBJEXT) \
  m_sax_dictionary.$(LTOBJEXT) \
  m_xml_error.$(LTOBJEXT) \
  m_xml_parser.$(LTOBJEXT)

