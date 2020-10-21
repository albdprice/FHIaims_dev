# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the CML examples of XMLF90
#

AM_CPPFLAGS = -I../../../src/cml -I../../../src/wxml
AM_FCFLAGS = -I../../../src/cml -I../../../src/wxml

noinst_PROGRAMS = example

example_SOURCES = example.f90
example_LDADD = ../../../src/libxmlf90.la

TESTS = example

CLEANFILES = *.mod output.xml
