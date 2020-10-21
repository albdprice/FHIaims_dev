# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the XPATH examples of XMLF90
#

AM_CPPFLAGS = -I../../../src/sax -I../../../src/xpath -DXPATH_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../src/sax -I../../../src/xpath

noinst_PROGRAMS = pseudo relative text

pseudo_SOURCES = pseudo.F90 m_pseudo_types.f90
pseudo_LDADD = ../../../src/libxmlf90.la
pseudo_DEPENDENCIES = pseudo.xml

relative_SOURCES = relative.F90
relative_LDADD = ../../../src/libxmlf90.la
relative_DEPENDENCIES = Ba.xml

text_SOURCES = text.F90
text_LDADD = ../../../src/libxmlf90.la
text_DEPENDENCIES = Ba.xml

if F90_MOD_UPPERCASE
f03_mods = M_PSEUDO_TYPES.$(MODEXT)
else
f03_mods = m_pseudo_types.$(MODEXT)
endif

EXTRA_DIST = Ba.xml README pseudo.xml

CLEANFILES = $(f03_mods)

TESTS = pseudo relative text

pseudo.$(OBJEXT): m_pseudo_types.$(OBJEXT)
