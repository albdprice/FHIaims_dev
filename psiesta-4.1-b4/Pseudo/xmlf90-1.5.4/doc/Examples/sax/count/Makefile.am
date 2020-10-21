# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the SAX examples of XMLF90
#

# Connect this directory to the main source code
AM_CPPFLAGS = -I../../../../src/sax -DSAX_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../../src/sax

noinst_PROGRAMS = count

count_SOURCES = count.F90 m_count.f90
count_LDADD = ../../../../src/libxmlf90.la
count_DEPENDENCIES = big-file.xml

if F90_MOD_UPPERCASE
f03_mods = M_COUNT.$(MODEXT)
else
f03_mods = m_count.$(MODEXT)
endif

EXTRA_DIST = README big-file.xml

CLEANFILES = $(f03_mods)

TESTS = count

count.$(OBJEXT): m_count.$(OBJEXT)
