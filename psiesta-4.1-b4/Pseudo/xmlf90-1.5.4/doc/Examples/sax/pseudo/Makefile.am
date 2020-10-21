# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the SAX examples of XMLF90
#

# Connect this directory to the main source code
AM_CPPFLAGS = -I../../../../src/sax -DSAX_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../../src/sax

noinst_PROGRAMS = pseudo

pseudo_SOURCES = pseudo.F90 m_pseudo.f90 m_pseudo_types.f90
pseudo_LDADD = ../../../../src/libxmlf90.la
pseudo_DEPENDENCIES = pseudo.xml

if F90_MOD_UPPERCASE
f03_mods = M_PSEUDO.$(MODEXT) M_PSEUDO_TYPES.$(MODEXT)
else
f03_mods = m_pseudo.$(MODEXT) m_pseudo_types.$(MODEXT)
endif

EXTRA_DIST = README pseudo.xml

CLEANFILES = $(f03_mods)

TESTS = pseudo

pseudo.$(OBJEXT): m_pseudo.$(OBJEXT)
m_pseudo.$(OBJEXT): m_pseudo_types.$(OBJEXT)
