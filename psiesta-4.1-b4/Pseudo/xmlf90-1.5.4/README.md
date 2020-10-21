XMLF90 README
=============

xmlf90 is a suite of libraries to handle XML in Fortran. It has two
major components:

- A XML parsing library. The parser was designed to be a useful
  tool in the extraction and analysis of data in the context of
  scientific computing, and thus the priorities were efficiency and the
  ability to deal with very large XML files while maintaining a small
  memory footprint. The most complete programming interface is
  based on the very successful SAX (Simple API for XML) model,
  although a partial DOM interface and a very experimental XPATH interface
  are also present.
- A library (xmlf90-wxml) that facilitates the writing of well-formed
  XML, including such features as automatic start-tag completion,
  attribute pretty-printing, and element indentation. There are also
  helper routines to handle the output of numerical arrays.

NOTE: The FoX project, started by Toby White and now maintained by
Andrew Walker, has produced a more robust and feature-rich package
starting from xmlf90, although not as fast and trim in some areas.  It
can be found in: http://www1.gly.bris.ac.uk/~walker/FoX/


INSTALLATION
------------

To install the library, just use the standard procedure:

    ./configure --prefix=/path/to/install/directory
    make
    make check
    make install

You can go into the subdirectories 'doc/Examples' and explore, and go
into 'doc/Tutorial' and try the exercises in the User Guide.


COMPILING USER PROGRAMS
-----------------------

After installation, the appropriate modules and library files should
already be in \$PREFIX/include and \$PREFIX/lib, respectively, where
\$PREFIX is the installation directory specified in the --prefix option
above.

To compile user programs, it is suggested that the user create a
separate directory to hold the program files and prepare a Makefile
following this example (FC, FFLAGS, and LDFLAGS need to be set appropriately
for the Fortran compiler used):

    #---------------------------------------------------------------
    #
    default: example
    #
    #---------------------------
    XMLF90_ROOT=/path/to/installation
    #
    XMLF90_LIBS=-L$(XMLF90_ROOT)/lib -lxmlf90
    XMLF90_INCFLAGS=-I$(XMLF90_ROOT)/include
    #
    INCFLAGS := $(XMLF90_INCFLAGS) $(INCFLAGS)
    LIBS:= $(XMLF90_LIBS) $(LIBS)
    #---------------------------
    #
    OBJS= m_handlers.o example.o
     
    example:  $(OBJS)
            $(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
    #
    clean: 
            rm -f *.o example *.mod
    #
    # Building rules
    #
    .F.o:
            $(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
    .f.o:
            $(FC) -c $(FFLAGS) $(INCFLAGS)   $<
    .F90.o:
            $(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
    .f90.o:
            $(FC) -c $(FFLAGS) $(INCFLAGS)   $<
    #
    #---------------------------------------------------------------

Here it is assumed that the user has two source files,
'example.f90' and 'm_handlers.f90'. Simply typing
'make' will compile 'example', pulling in all the needed
modules and library objects.

Alternatively, the two lines defining explicitly the XMLF90_LIBS
and XMLF90_INCFLAGS variables can be replaced by the line:

include \$(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk

which will perform that function abstracting the details from the user.

As of version 1.5.3 of xmlf90, the xmlf90.mk will force the use of
static libraries, so that no other symbols but XMLF90_ROOT need to be
exported for programs to work.

