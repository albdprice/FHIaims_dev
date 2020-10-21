# -*- Autoconf -*-
#
# Configure script for the XMLF90 package
#
# Copyright (C) 2017 Y. Pouillon
#
# This file is part of the XMLF90 software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Autotools startup
#

# Init Autoconf
AC_PREREQ(2.69)
AC_INIT([XMLF90], [1.5.4],
  [https://bugs.launchpad.net/xmlf90],
  [xmlf90],
  [https://launchpad.net/xmlf90])
AC_REVISION([Autotools support for XMLF90])
AC_CONFIG_AUX_DIR(config/gnu)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR(src/sax/m_xml_parser.f90)

# Make relevant path available globally
_AC_SRCDIRS([.])
AC_SUBST(ac_abs_top_srcdir)
AC_SUBST(ac_abs_top_builddir)

# Init Automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.14])

# Generate timestamp
x90_timestamp=`date '+%Y%m%dT%H%M%S%z'`
AC_SUBST(x90_timestamp)

# Export version triplet
XMLF90_VERSION_MAJOR=`echo "${PACKAGE_VERSION}" | cut -d. -f1`
XMLF90_VERSION_MINOR=`echo "${PACKAGE_VERSION}" | cut -d. -f2`
XMLF90_VERSION_MICRO=`echo "${PACKAGE_VERSION}" | cut -d. -f3`
AC_SUBST(XMLF90_VERSION_MAJOR)
AC_SUBST(XMLF90_VERSION_MINOR)
AC_SUBST(XMLF90_VERSION_MICRO)

                    # ------------------------------------ #

#
# Command-line options
#

# Control the generation of documentation
AC_ARG_ENABLE([docgen],
  AC_HELP_STRING([--enable-docgen],
    [Enable the generation of developer documentation (experimental, disabled by default)]),
  [
    if test "${enableval}" = "no" -o "${enableval}" = "yes"; then
      x90_docgen_enable="${enableval}"
    else
      AC_MSG_ERROR([--enable-docgen must be 'yes' or 'no', not '${enableval}'])
    fi
  ],
  [x90_docgen_enable="no"])
AC_SUBST(x90_docgen_enable)

# Allow the user to specify where to find Doxygen
AC_ARG_VAR([DOXYGEN], [Path to the 'doxygen' command])

                    # ------------------------------------ #

#
# System utilities
#

# Check for common programs
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP

                    # ------------------------------------ #

#
# C language support
#

# Init C compiler and preprocessor
AC_PROG_CC
AC_PROG_CPP

# C compiler peculiarities (for Libtool)
AM_PROG_CC_C_O

                    # ------------------------------------ #

#
# Fortran language support
#

# Look for the Fortran compiler
if test "${FC}" != "" -a ! -x "${FC}"; then
  x90_fc_probe=`echo "${FC}" | sed -e 's/ .*//'`
  if test ! -x "${x90_fc_probe}"; then
    AC_PATH_PROG([x90_fc_path],[${x90_fc_probe}])
    if test "${x90_fc_path}" = ""; then
      AC_MSG_ERROR([could not run Fortran compiler "${FC}"])
    fi
  fi
fi
AC_PROG_FC

# Fail if no Fortran compiler is available
if test "${FC}" = ""; then
  AC_MSG_ERROR([no Fortran compiler available])
fi

# Fortran compiler peculiarities
AX_F90_MODULE_EXTENSION
AX_F90_MODULE_CASE
MODEXT="${ax_cv_f90_modext}"
AC_SUBST(MODEXT)

# Select how to handle Fortran modules
AM_CONDITIONAL([F90_MOD_UPPERCASE],
  [test "${ax_cv_f90_mod_uppercase}" = "yes"])

                    # ------------------------------------ #

#
# Libtool configuration
#

# Init Libtool (must be done once compilers are fully set)
LT_INIT
LT_PREREQ([2.4.2])
LT_LANG([Fortran])
LTOBJEXT="lo"
AC_SUBST(LTOBJEXT)

                    # ------------------------------------ #

#
# Documentation generation
#

# Look for Doxygen
if test "${x90_docgen_enable}" != "no"; then
  if test -z "${DOXYGEN}"; then
    AC_CHECK_PROGS([DOXYGEN], [doxygen])
  fi
fi

# Look for LaTeX
#
# Note: These environment variables must be reset because Doxygen uses
#       hard-coded program names.
#
MAKEINDEX=""
PDFLATEX=""
if test "${x90_docgen_enable}" != "no" -a ! -z "${DOXYGEN}"; then
  AC_CHECK_PROGS([MAKEINDEX], [makeindex])
  AC_CHECK_PROGS([PDFLATEX], [pdflatex])
fi
AC_SUBST(MAKEINDEX)
AC_SUBST(PDFLATEX)

                    # ------------------------------------ #

#
# Final checks
#

# Write down YAML configuration now (allows full report with deferred errors)
AC_OUTPUT([config/data/xmlf90-config.yml])

# Decide what to do with documentation
if test "${x90_docgen_enable}" = "yes" -a -z "${DOXYGEN}"; then
  AC_MSG_FAILURE([Documentation generation was requested from the command line
                  but Doxygen was not found])
fi
if test "${x90_docgen_enable}" = "auto" -a ! -z "${DOXYGEN}"; then
  x90_docgen_enable="yes"
fi
AM_CONDITIONAL([DO_BUILD_DOXYGEN], [test "${x90_docgen_enable}" = "yes"])
AM_CONDITIONAL([DO_BUILD_LATEX], [test "${x90_docgen_enable}" = "yes" -a ! -z "${MAKEINDEX}" -a ! -z "${PDFLATEX}"])

                    # ------------------------------------ #

#
# Output configuration
#

AC_OUTPUT([
  config/data/xmlf90.mk
  config/data/xmlf90.pc
  Makefile
  src/Makefile
  src/xmlf90_info.f90
  src/cml/Makefile
  src/dom/Makefile
  src/sax/Makefile
  src/wxml/Makefile
  src/xpath/Makefile
  doc/Doxyfile
  doc/Makefile
  doc/Examples/Makefile
  doc/Examples/cml/Makefile
  doc/Examples/dom/Makefile
  doc/Examples/sax/Makefile
  doc/Examples/sax/count/Makefile
  doc/Examples/sax/features/Makefile
  doc/Examples/sax/pseudo/Makefile
  doc/Examples/sax/simple/Makefile
  doc/Examples/wxml/Makefile
  doc/Examples/xpath/Makefile
  doc/Tutorial/Makefile
  doc/Tutorial/sax/Makefile
  doc/Tutorial/xpath/Makefile
])
                    # ------------------------------------ #

#
# Report configuration
#

AC_MSG_NOTICE([

Final build parameters
----------------------

    TSTAMP     = ${x90_timestamp}

    docgen     = ${x90_docgen_enable}

    CPP        = ${CPP}
    CPPFLAGS   = ${CPPFLAGS}
    CC         = ${CC}
    CFLAGS     = ${CFLAGS}
    FC         = ${FC}
    FCFLAGS    = ${FCFLAGS}
    LDFLAGS    = ${LDFLAGS}
    LIBS       = ${LIBS}

    DOXYGEN    = ${DOXYGEN}
    MAKEINDEX  = ${MAKEINDEX}
    PDFLATEX   = ${PDFLATEX}

This information is also available in YAML format inside
config/data/${PACKAGE}-config.yml for future reference.

You can now build ${PACKAGE_STRING} by typing "make", or
e.g. "make -j8" for a parallel build with 8 processes.
])
