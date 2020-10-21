# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the SAX examples of XMLF90
#

# Connect this directory to the main source code
AM_CPPFLAGS = -I../../../../src/sax -DSAX_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../../src/sax

noinst_PROGRAMS = example xmlcheck

example_SOURCES = example.F90 m_handlers.f90
example_LDADD = ../../../../src/libxmlf90.la
example_DEPENDENCIES = test.xml

xmlcheck_SOURCES = xmlcheck.f90
xmlcheck_LDADD = ../../../../src/libxmlf90.la

if F90_MOD_UPPERCASE
f03_mods = M_HANDLERS.$(MODEXT)
else
f03_mods = m_handlers.$(MODEXT)
endif

EXTRA_DIST = README test.xml

CLEANFILES = $(f03_mods)

TESTS = example xmlcheck

example.$(OBJEXT): m_handlers.$(OBJEXT)
