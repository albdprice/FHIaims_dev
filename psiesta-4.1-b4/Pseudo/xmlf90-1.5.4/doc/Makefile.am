# Copyright (C) 2017 Yann Pouillon

#
# Makefile for the documentation of XMLF90
#

SUBDIRS = Examples Tutorial

EXTRA_DIST = \
  ReleaseNotes-1.1 \
  ReleaseNotes-1.2 \
  ReleaseNotes-1.4 \
  ReleaseNotes-1.5

if DO_BUILD_DOXYGEN
if DO_BUILD_LATEX
all-local: .docgen-stamp .docgen-latex-stamp
else
all-local: .docgen-stamp
endif

.docgen-stamp: Doxyfile
	$(DOXYGEN) Doxyfile
	touch .docgen-stamp
if DO_BUILD_LATEX
developers/@PACKAGE@-@VERSION@-reference.pdf: .docgen-latex-stamp

.docgen-latex-stamp: .docgen-stamp
	cd developers/latex && $(MAKE) all
	$(LN_S) latex/refman.pdf developers/@PACKAGE@-@VERSION@-reference.pdf
	touch .docgen-latex-stamp
endif

install-data-local:
	$(INSTALL) -m 755 -d $(DESTDIR)$(docdir)/developers
	find developers/html -type f -exec $(INSTALL_DATA) -D "{}" "$(DESTDIR)$(docdir)/{}" \;
if DO_BUILD_LATEX
	$(INSTALL_DATA) developers/@PACKAGE@-@VERSION@-reference.pdf $(DESTDIR)$(docdir)/developers
endif

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)/developers

clean-local:
	rm -rf developers
	rm -f .docgen-stamp .docgen-latex-stamp
endif
