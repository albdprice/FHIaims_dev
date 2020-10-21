# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the DOM examples of XMLF90
#

# Connect this directory to the main source code
AM_CPPFLAGS = -I../../../src/dom -I../../../src/sax -DDOM_DATA_DIR=\"$(srcdir)/\"
AM_FCFLAGS = -I../../../src/dom -I../../../src/sax

# Example programs
noinst_PROGRAMS = data features pseudo pseudo_dom text

data_SOURCES = data.F90
data_LDADD = ../../../src/libxmlf90.la
data_DEPENDENCIES = pseudo.xml

features_SOURCES = features.F90
features_LDADD = ../../../src/libxmlf90.la
features_DEPENDENCIES = test.xml

pseudo_SOURCES = pseudo.F90 m_pseudo_types.f90
pseudo_LDADD = ../../../src/libxmlf90.la
pseudo_DEPENDENCIES = pseudo.xml

pseudo_dom_SOURCES = pseudo_dom.F90 m_pseudo_types.f90 m_psdom.f90
pseudo_dom_LDADD = ../../../src/libxmlf90.la
pseudo_dom_DEPENDENCIES = pseudo.xml

text_SOURCES = text.F90
text_LDADD = ../../../src/libxmlf90.la
text_DEPENDENCIES = big-file.xml

# Fortran modules
if F90_MOD_UPPERCASE
f03_mods = M_PSDOM.$(MODEXT) M_PSEUDO_TYPES.$(MODEXT)
else
f03_mods = m_psdom.$(MODEXT) m_pseudo_types.$(MODEXT)
endif

# Data required to run the programs
xml_data = \
  big-file.xml \
  pseudo.xml \
  test.xml

# Make sure all important files are distributed
EXTRA_DIST = \
  README \
  $(xml_data)

TESTS = data features pseudo pseudo_dom text

# Fortran modules must be explicitly removed
CLEANFILES = \
  $(f03_mods) \
  all.xml \
  clone.xml \
  dumptext.xml \
  features.xml \
  node.xml \
  nosemilocal.xml \
  out.xml \
  ref.xml

# Explicit dependencies
m_psdom.$(OBJEXT): m_pseudo_types.$(OBJEXT)
pseudo.$(OBJEXT): m_pseudo_types.$(OBJEXT)
pseudo_dom.$(OBJEXT): m_pseudo_types.$(OBJEXT) m_psdom.$(OBJEXT)
