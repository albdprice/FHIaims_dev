Fortran XML Tools
=================

SAX
---

Flib SAX is a SAX level 1.0 implementation in Fortran 90.

A PDF Tutorial and UserGuide is available [here](UserGuide.pdf)


Stream Xpath
------------

Stream Xpath is a library that emulates some of the features of the
Xpath standard, but working within the stream model of SAX.

Its small memory footprint makes it quite useful to process large
datafiles, for which the standard Xpath (built on top of the
memory-intensive DOM) would not be appropriate. However, the stream
paradigm forces the user to be careful about controlling the state of
the parser.

A PDF Tutorial and UserGuide is available [here](UserGuide.pdf)


WXML
----

WXML is a library that facilitates the writing of well-formed XML,
including such features as automatic start-tag completion, attribute
pretty-printing, and element indentation. There are also helper routines
to handle the output of numerical arrays.

[Desription of the routines](WXML.html)

See also the examples in the `Examples/wxml` subdirectory of the main
distribution.

Jon Wakelin has written a CML-formatting library on top of a slightly
modified WXML. Documentation is available [here](jumbo.html). For
examples of CML-formatting in strict WXML, see the `Examples/cml`
subdirectory of the main xmlf90 distribution. The two strands of WXML
will be merged very soon.


FDOM
----

FDOM is a a DOM level 1.0 implementation in Fortran 95. We have
implemented almost all the instance methods, although it is unlikely
that any of the class methods will ever be implemented. The FDOM is
still evolving but is already in a usable state. More importantly, as
all of the interfaces are standard, changes to the code will only take
place behind the scenes.

A page containing a breakdown of the FDOM methods is available
[here](DOM.html)

See also the examples in the `Examples/dom` subdirectory of the main
distribution.

\

*Jon Wakelin, Alberto Garcia, April 2004*
------------------------------------------------------------------------


Guidelines for developers
-------------------------

The parser is built on several levels:

1. Upper-level modules
  * m\_xml\_parser: The main module
  * m\_error     : Basic error handling
2. Intermediate layer
  * m\_sax\_fsm (A finite-state machine to parse the input)
3. Basic data structures and file interfaces
  * m\_sax\_reader:  File interface and character handling as per XML specs.
  * m\_sax\_buffer:  Basic homemade "variable length string", with some
  limitations (size, of course), but avoiding the use of dynamic structures
  for now.
  * m\_sax\_dictionary: Simple, not dynamic.
  * m\_sax\_charset: A simple hashing method for sets of characters.
  * m\_sax\_elstack: Simple stack to check well-formedness.
  * m\_sax\_entities: Entity replacement utilities.
4. Something which does not really belong in the parser but which
  is useful to massage the data extracted from the file:
  * m\_sax\_converters: Routines to turn pcdata chunks into numerical arrays

There are obviously a number of hardwired limitations, which should be
removed in a later version:

* Buffer size in buffer\_t definition. This is not as serious as it
  looks. Only long unbroken pcdata sections and overly long attribute
  names or values will be affected. Long SGML declarations and comments
  might be truncated, but they are not relevant anyway.

* Maximum number of attributes in an element tag (set in m\_sax\_dictionary)

While the parser does not use any variable-length strings (to keep it
compatible with existing Fortran90 compilers) or dynamical data
structures for attribute dictionaries, etc, such improvements could be
incorporated almost as drop-in replacements for existing sub-modules.

The coding style is that of the F subset of Fortran90. I strongly
believe that it makes for better coding and fewer errors.
Go to http:\/\/www.fortran.com/imagine1/ and get a feel for it. You can
download free implementations for Linux and Windows, or get an
inexpensive CD+Book combination to help support the project. Of course,
F *is* Fortran, so you can always compile it with a Fortran compiler.

