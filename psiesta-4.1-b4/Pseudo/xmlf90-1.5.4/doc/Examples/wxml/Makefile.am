# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the WXML examples of XMLF90
#

AM_CPPFLAGS = -I../../../src/wxml
AM_FCFLAGS = -I../../../src/wxml

noinst_PROGRAMS = pseudo simple

pseudo_SOURCES = pseudo.f90 m_pseudo_utils.f90
pseudo_LDADD = ../../../src/libxmlf90.la
pseudo_DEPENDENCIES = C.psf Fe.psf

simple_SOURCES = simple.f90
simple_LDADD = ../../../src/libxmlf90.la

if F90_MOD_UPPERCASE
f03_mods = M_PSEUDO_UTILS.$(MODEXT)
else
f03_mods = m_pseudo_utils.$(MODEXT)
endif

EXTRA_DIST = C.psf Fe.psf

CLEANFILES = $(f03_mods) C.xml Fe.xml simple.xml

TESTS = simple

check-local:
	echo "Testing pseudo with C.psf"
	cp $(srcdir)/C.psf PSF
	chmod u+w PSF
	./pseudo
	test -s PSXML
	rm PSF
	mv PSXML C.xml
	echo "Testing pseudo with Fe.psf"
	cp $(srcdir)/Fe.psf PSF
	chmod u+w PSF
	./pseudo
	test -s PSXML
	rm PSF
	mv PSXML Fe.xml

pseudo.$(OBJEXT): m_pseudo_utils.$(OBJEXT)
