# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the SAX tutorials of XMLF90
#

AM_CPPFLAGS = -I../../../src/sax -DSAX_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../src/sax

noinst_PROGRAMS = inventory simple

inventory_SOURCES = inventory.F90 m_handlers.f90
inventory_LDADD = ../../../src/libxmlf90.la
inventory_DEPENDENCIES = inventory.xml

simple_SOURCES = simple.F90
simple_LDADD = ../../../src/libxmlf90.la
simple_DEPENDENCIES = inventory.xml

if F90_MOD_UPPERCASE
f03_mods = M_AUX.$(MODEXT) M_HANDLERS.$(MODEXT)
else
f03_mods = m_aux.$(MODEXT) m_handlers.$(MODEXT)
endif

CLEANFILES = $(f03_mods)

EXTRA_DIST = README inventory.xml

TESTS = inventory simple

inventory.$(OBJEXT): m_handlers.$(OBJEXT)
