# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the XPATH tutorials of XMLF90
#

AM_CPPFLAGS = -I../../../src/sax -I../../../src/xpath -DXPATH_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../src/sax -I../../../src/xpath

noinst_PROGRAMS = \
  context \
  euros \
  item_context \
  item_context2 \
  mixing \
  self \
  simple \
  twoelements

context_SOURCES = context.F90
context_LDADD = ../../../src/libxmlf90.la
context_DEPENDENCIES = inventory.xml

euros_SOURCES = euros.F90
euros_LDADD = ../../../src/libxmlf90.la
euros_DEPENDENCIES = inventory.xml

item_context_SOURCES = item_context.F90
item_context_LDADD = ../../../src/libxmlf90.la
item_context_DEPENDENCIES = inventory.xml

item_context2_SOURCES = item_context2.F90
item_context2_LDADD = ../../../src/libxmlf90.la
item_context2_DEPENDENCIES = inventory.xml

mixing_SOURCES = mixing.F90
mixing_LDADD = ../../../src/libxmlf90.la
mixing_DEPENDENCIES = inventory_text.xml

self_SOURCES = self.F90
self_LDADD = ../../../src/libxmlf90.la
self_DEPENDENCIES = inventory.xml

simple_SOURCES = simple.F90
simple_LDADD = ../../../src/libxmlf90.la
simple_DEPENDENCIES = inventory.xml

twoelements_SOURCES = twoelements.F90
twoelements_LDADD = ../../../src/libxmlf90.la
twoelements_DEPENDENCIES = inventory.xml

if F90_MOD_UPPERCASE
f03_mods = M_AUX.$(MODEXT)
else
f03_mods = m_aux.$(MODEXT)
endif

CLEANFILES = $(f03_mods)

EXTRA_DIST = README inventory.xml inventory_text.xml

TESTS = context euros item_context item_context2 mixing self simple twoelements
