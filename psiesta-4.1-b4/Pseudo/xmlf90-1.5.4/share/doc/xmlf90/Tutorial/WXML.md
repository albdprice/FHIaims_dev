WXML Library
============

## Routines

### General routines

1.  [xml\_OpenFile](#init) - Mandatory Initialization routine
2.  [xml\_Close](#close) - Mandatory finalization routine, closes
    channels, etc
3.  [str](#str) - utility to convert reals and integers to character
    strings

### XML routines

1.  [xml\_NewElement](#start) - writes an xml start tag
2.  [xml\_AddAttribute](#attr) - adds an attribute to a tag
3.  [xml\_AddPcdata](#text) - adds text to an xml element
4.  [xml\_AddArray](#array) - dumps the contents of an array as pcdata
5.  [xml\_EndElement](#end) - writes an xml end tag

Subroutine Guide

1.  #### [xml\_OpenFile(filename, ind, xf)]()

      argument   role                             type                optional   default
      ---------- -------------------------------- ------------------- ---------- ---------
      filename   xml filename                     character(len=\*)   no          
      ind        controls indentation of output   logical             yes        .true.
      xf         xml filename                     type(xmlf\_t)       no          

2.  #### [xml\_Close(filehandle)]()

      argument   role             type            optional   default
      ---------- ---------------- --------------- ---------- ---------
      xf         xml filehandle   type(xmlf\_t)   no          

3.  #### [function str(value, format)]()

      argument   role                         type                                 optional   default
      ---------- ---------------------------- ------------------------------------ ---------- ---------
      value      value to convert to string   real\*8, real\*4, integer, logical   no          
      format     format for reals             character(len=\*)                    yes        g22.12

    ------------------------------------------------------------------------

<!-- -->

1.  #### [xml\_NewElement(xf, name)]()

      argument   role                 type                optional   default
      ---------- -------------------- ------------------- ---------- ---------
      xf         xml filehandle       type(xmlf\_t)       no          
      name       name of tag to add   character(len=\*)   no          

2.  #### [xml\_AddAttribute(xf, attname, value)]()

      argument   role              type                                      optional   default
      ---------- ----------------- ----------------------------------------- ---------- ---------
      xf         xml filehandle    type(xmlf\_t)                             no          
      attname    attribute name    character(len=\*)                         no          
      value      attribute value   character(len=\*) (convert using str())   no          

3.  #### [xml\_AddPcdata(xf, pcdata)]()

      argument   role             type                                              optional   default
      ---------- ---------------- ------------------------------------------------- ---------- ---------
      xf         xml filehandle   type(xmlf\_t)                                     no          
      pcdata     string to add    character(len=\*) (convert numbers using str())   no          

4.  #### [xml\_AddArray(xf, a, format)]()

      argument   role             type                    optional   default
      ---------- ---------------- ----------------------- ---------- ---------------------
      xf         xml filehandle   type(xmlf\_t)           no          
      a          array (:)        integer, real, double   no          
      format     format           character(len=\*)       yes        6(i12) / 4(es20.12)

5.  #### [xml\_EndElement(xf, name)]()

      argument   role                       type                optional   default
      ---------- -------------------------- ------------------- ---------- ---------
      xf         xml filehandle             type(xmlf\_t)       no          
      name       name of element to close   character(len=\*)   no          

------------------------------------------------------------------------

## Jumbo90

Jumbo90 is a library written in Fortran for creating CML documents.
Actually, Jumbo allows you to create XML, CML and STMML (another markup
language closely related to CML). After a moments thought you may ask
why do we need a library to write things? I mean, you can always type

    write(*,*)
     "<cml>stuff</cml>"

Well, while in theory it is a fairly simply task to write CML, or indeed
any XML, in any language, in practice it's not quite that easy.
Balancing tags that are opened at one end of a 2000 line source file and
closed at the other or even opened and closed in different source files,
is an unnecessary memory exercise, and balancing them in the right order
a pain. More importantly, you'll probably only spot your mistake when
you try to process the file, and given that XML is case sensitive it's a
pretty easy to make a mistake. Furthermore, CML elements will often be
built up in a standard manner using a large number of tags - and it can
be fairly tedious typing in all of these. Add to this the fact that it
is extremely difficult to correctly indent an XML document whose
structure is as dynamic as that of the underlying program \[1\] (just
using write statements). So, we decided that it was probably worth our
while to create a library for handling the more tedious or tricky
elements of creating a CML document. We also added functionality to
check the well-formedness of all XML. In fact, most programming
languages provide utilties for performing these seemingly simple tasks

The original idea for Jumbo90 came from Peter Murray-Rust, who has a
similar but much more rounded suite of tools written in Java, called
Jumbo. I then converted Peter's original F77 code into F90, but it
suffered a little from its F77 inheritance, so now the base XML
formatting layer of Jumbo is provided by Alberto Garcia's WXML library.

The Jumbo90 README containing a breakdown of all Jumbo90 subroutines is
available [here](jumbodocs.html)

\

------------------------------------------------------------------------

*Jon Wakelin, March 2004*

\[1\] It's only a cosmetic feature but it's very useful when you need
someone to read the file, of course, there's nothing wrong with writing
a 12M XML file in a single line.


### What is Jumbo90?

-   Jumbo90 is a CML formatting Library for Fortran, it provides:
-   Convenience routines to write complete CML elements
-   Convenience routines to write complete STMML elements
-   Optional Indentation of output
-   Checks that tags are properly balanced
-   Checks that attribute names are well formed (i.e. contain only \[A-Z
    a-z 0-9 \_ -\] and start with \[A-Z a-z \_\]).
    -   The colon is allowed but should be reserved for namespacing so
        an error is flagged if more than one colon is present in
        attribute name
-   Checks if attribute values contain predefined entities (&lt;, &gt;,
    &, ', ")
-   Checks for duplicate attribute names

How to use Jumbo90

Jumbo covers CML, STMML and basic XML. Below is the full list of
subroutine names, followed by the arguments they take. Whenever writing
a real\*8 or real\*4 you can always pass an optional format argument,
e.g.

        
        call xml_AddAttribute(file, 'martin','height', 32.1285)
        call xml_AddAttribute(file, 'martin','height', 32.1285, '(f6.2)')

would add an "height" attribute with a value of "32.1285" to a "martin"
element . In the first case the default format '(f8.3)' would be used,
in the second case the user supplied format '(f6.2)' would be used. Many
subroutines can take a number of optional arguments (reflecting the CML
schema) therefore the longer 'argument=value' format should be use when
calling these subroutines, e.g.,


        call xml_OpenFile('myFile.xml', file)

Notes:
------

Jumbo90 is based on an earlier fortran 77 program called Jumbo77, which
was itself based on an existing Java CML parser called JUMBO written by
*Peter Murray-Rust*. Jumbo90 is now built in a modular fashion allowing
output of basic XML, STMML and CML. The STMML layer builds on the XML.
The CML layer, in turn, builds on the STMML layer. The base XML writing
utilities are now provided by (xmlf90-wxml), a set of F90 modules
written by *Alberto Garcia*.


### General routines

1.  [xml\_OpenFile](#init) - Mandatory Initialization routine
2.  [xml\_Close](#close) - Mandatory finalization routine, closes
    channels, etc
3.  [str](#str) - utility to convert, floats and integers to character
    strings

### XML routines

1.  [xml\_NewElement](#start) - writes an xml start tag
2.  [xml\_AddAttribute](#attr) - adds an attribute to a tag
3.  [xml\_AddPcdata](#text) - adds text to an xml element
4.  [xml\_EndElement](#end) - writes an xml end tag

### STMML routines

1.  [stmAddStartTag](#stmstart)
2.  [stmAddScalar](#stmsca)
3.  [stmAddArray](#stmarr)
4.  [stmAddMatrix](#stmmat)
5.  [stmAddTriangle](#stmtri)
6.  [stmError](#stmerr)
7.  [stmMessage](#stmmsg)
8.  [stmWarning](#stmwar)

### CMLCore routines

1.  [cmlAddMolecule](#cmlmol) - adds a complete CML &lt;molecule&gt;
    element
2.  [cmlAddAtom](#cmlato) - adds a CML &lt;atom&gt; start tag
3.  [cmlAddCoordinates](#cmlcoo) - adds coordinate attributes to an
    &lt;atom&gt; tag
4.  [cmlAddCrystal](#cmlcry) - adds a complete CML &lt;crystal&gt;
    element
5.  [cmlAddMetadata](#cmlmet) - adds a complete CML &lt;metadata&gt;
    element
6.  [cmlAddLength](#cmllen)(length, id, atomRef1, atomRef2, fmt)
7.  [cmlAddAngle](#cmlang)(angle, id, atomRef1, atomRef2, atomRef3, fmt)
8.  [cmlAddTorsion](#cmltor)(torsion, id, atomRef1, atomRef2, atomRef3,
    atomRef4, fmt)
9.  [cmlAddEigenvalue](#cmleig)(n, dim, eigvec, eigval, id, title,
    dictRef, fmt)

### CMLComa routines (Condensed Matter)

1.  [cmlAddLattice](#comalat) - write a complete lattice element
2.  [cmlAddProperty](#comapro) - write a complete property element,
    containing scalar, array (vector) or matrix value (i.e. it has
    several interfaces)
3.  [cmlAddParameter](#comapar) - adds a complete CML &lt;parameter&gt;
    element

Subroutine Guide

1.  #### [xml\_OpenFile(filename, ind, xf)]()

      argument   role                             type                optional   default
      ---------- -------------------------------- ------------------- ---------- ---------
      filename   xml filename                     character(len=\*)   no          
      ind        controls indentation of output   logical             yes        .true.
      xf         xml filename                     type(xmlf\_t)       no          

2.  #### [xml\_Close(filehandle)]()

      argument   role             type            optional   default
      ---------- ---------------- --------------- ---------- ---------
      xf         xml filehandle   type(xmlf\_t)   no          

3.  #### [function str(value, format)]()

      argument   role                         type                        optional   default
      ---------- ---------------------------- --------------------------- ---------- ---------
      value      value to convert to string   real\*8, real\*4, integer   no          
      format     format for reals             character(len=\*)           yes        g22.12

    ------------------------------------------------------------------------

<!-- -->

1.  #### [xml\_NewElement(xf, name)]()

      argument   role                 type                optional   default
      ---------- -------------------- ------------------- ---------- ---------
      xf         xml filehandle       type(xmlf\_t)       no          
      name       name of tag to add   character(len=\*)   no          

2.  #### [xml\_AddAttribute(xf, name, attname, value)]()

      argument   role               type                                              optional   default
      ---------- ------------------ ------------------------------------------------- ---------- ---------
      xf         xml filehandle     type(xmlf\_t)                                     no          
      name       name of tag        character(len=\*)                                 no          
      attname    attribute name     character(len=\*)                                 no          
      value      attribute value    character(len=\*) | integer | real\*8 | real\*4   no          
      fmt        format for reals   character(len=\*)                                 yes        f8.3

3.  #### [xml\_AddPcdata(xf, value)]()

      argument   role               type                                              optional   default
      ---------- ------------------ ------------------------------------------------- ---------- ---------
      xf         xml filehandle     type(xmlf\_t)                                     no          
      value      string to add      character(len=\*) | integer | real\*8 | real\*4   no          
      fmt        format for reals   character(len=\*)                                 yes        f8.3

4.  #### [xml\_EndElement(xf, name)]()

      argument   role                       type                optional   default
      ---------- -------------------------- ------------------- ---------- ---------
      xf         xml filehandle             type(xmlf\_t)       no          
      name       name of element to close   character(len=\*)   no          

------------------------------------------------------------------------

1.  #### [stmAddStartTag(xf, name, id, title, dictref)]()

      argument   role                   type                optional   default
      ---------- ---------------------- ------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)       no          
      name       tag name               character(len=\*)   no          
      id         unique id              character(len=\*)   yes         
      title      tag description        character(len=\*)   yes         
      dictref    dictionary reference   character(len=\*)   yes         

2.  #### [stmAddScalar(xf, value, id, title, dictref, type, fmt)]()

      argument   role                   type                                              optional   default
      ---------- ---------------------- ------------------------------------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)                                     no          
      value      the scalar value       character(len=\*) | integer | real\*8 | real\*4   no          
      id         unique id              character(len=\*)                                 yes         
      title      tag description        character(len=\*)                                 yes         
      dictref    dictionary reference   character(len=\*)                                 yes         
      fmt        format for reals       character(len=\*)                                 yes        f8.3

3.  #### [stmAddArray(xf, nvalue, array, id, title, dictref, type|units, delim|fmt)]()

      argument   role                   type                                              optional   default
      ---------- ---------------------- ------------------------------------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)                                     no          
      nvalue     length of array        integer                                           no          
      array      the array              character(len=\*) | integer | real\*8 | real\*4   no          
      id         unique id              character(len=\*)                                 yes         
      title      tag description        character(len=\*)                                 yes         
      dictref    dictionary reference   character(len=\*)                                 yes         
      fmt        format for reals       character(len=\*)                                 yes        f8.3

4.  #### [stmAddMatrix(xf, nrows, ncols, dim, matrix, id, title, dictref, units, fmt)]()

      argument   role                   type                                              optional   default
      ---------- ---------------------- ------------------------------------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)                                     no          
      nrows      number of rows         integer                                           no          
      ncols      number of columns      integer                                           no          
      dim        fastest dimension      integer                                           no          
      matrix     the matrix             character(len=\*) | integer | real\*8 | real\*4   no          
      id         unique id              character(len=\*)                                 yes         
      title      tag description        character(len=\*)                                 yes         
      dictref    dictionary reference   character(len=\*)                                 yes         
      fmt        format for reals       character(len=\*)                                 yes        f8.3

5.  #### [stmAddTriangle(xf, nvalue, array, id, title, dictref, units, fmt)]()

      argument   role                   type                                              optional   default
      ---------- ---------------------- ------------------------------------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)                                     no          
      nvalue     length of array        integer                                           no          
      array      the array              character(len=\*) | integer | real\*8 | real\*4   no          
      id         unique id              character(len=\*)                                 yes         
      title      tag description        character(len=\*)                                 yes         
      dictref    dictionary reference   character(len=\*)                                 yes         
      fmt        format for reals       character(len=\*)                                 yes        f8.3

6.  #### [stmAddError(xf, msg, id, title, dictRef)]()

      argument   role                   type                optional   default
      ---------- ---------------------- ------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)       no          
      msg        length of array        character(len=\*)   no          
      id         unique id              character(len=\*)   yes         
      title      tag description        character(len=\*)   yes         
      dictref    dictionary reference   character(len=\*)   yes         

7.  #### [stmAddMessage(xf, msg, id, title, dictRef)]()

      argument   role                   type                optional   default
      ---------- ---------------------- ------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)       no          
      msg        length of array        character(len=\*)   no          
      id         unique id              character(len=\*)   yes         
      title      tag description        character(len=\*)   yes         
      dictref    dictionary reference   character(len=\*)   yes         

8.  #### [stmAddWarning(xf, msg, id, title, dictRef)]()

      argument   role                   type                optional   default
      ---------- ---------------------- ------------------- ---------- ---------
      xf         xml filehandle         type(xmlf\_t)       no          
      msg        length of array        character(len=\*)   no          
      id         unique id              character(len=\*)   yes         
      title      tag description        character(len=\*)   yes         
      dictref    dictionary reference   character(len=\*)   yes         

------------------------------------------------------------------------

#### [cmlAddMolecule(xf, natoms, elements, coords, style, id, title, dictref, fmt)]()

  argument            role                                                       type                optional   default
  ------------------- ---------------------------------------------------------- ------------------- ---------- ---------
  xf                  xml filehandle                                             type(xmlf\_t)       no          
  natoms              number of atoms                                            integer             no          
  elements(natoms)    list of atomic symbols                                     character(len=2)    no          
  coords(3, natoms)   atomic coordinates                                         real\*8 | real\*4   no          
  style               [CML output style](#style) (x3 | xFrac | xyz3 | xyzFrac)   character(len=\*)   yes        x3
  id                  unique id                                                  character(len=\*)   yes         
  title               tag description                                            character(len=\*)   yes         
  dictref             dictionary reference                                       character(len=\*)   yes         
  fmt                 format for reals                                           character(len=\*)   yes        f8.3

#### [cmlAddAtom(xf, elem, id, charge, hCount, occupancy, fmt)]()

  argument    role                   type                optional   default
  ----------- ---------------------- ------------------- ---------- ---------
  xf          xml filehandle         type(xmlf\_t)       no          
  elem        atomic symbol          character(len=2)    yes         
  id          unique id              character(len=\*)   yes         
  charge      formal charge          integer             no          
  hCount      hydrogen count         integer             no          
  occupancy   site occupancy         real\*8 | real\*4   no          
  fmt         format for occupancy   character(len=\*)   yes        f8.3

#### [cmlAddCoordinates(xf, x, y, z, style, fmt)]()

  -----------------------------------------------------------------------------------------
  argument   role                                  type                optional   default
  ---------- ------------------------------------- ------------------- ---------- ---------
  xf         xml filehandle                        type(xmlf\_t)       no          

  x          *x* coordinate in cartesian format    real\*8 | real\*4   no          

  y          *y* coordinate in cartesian format    real\*8 | real\*4   no          

  z          *z* coordinate in cartesian format    real\*8 | real\*4   yes         

  style      [CML output style](#style)\           character(len=\*)   yes        x3
             (x3 | xFrac | xyz3 | xyzFrac | xy2)                                  

  fmt        format for coordinates                character(len=\*)   yes        f8.3
  -----------------------------------------------------------------------------------------

#### [cmlAddcrystal(xf, a, b, c, alpha, beta, gamma, id, title, dictref, lenunits, angunits, fmt)]()

  argument   role                           type                optional   default
  ---------- ------------------------------ ------------------- ---------- ----------
  xf         xml filehandle                 type(xmlf\_t)       no          
  a          lattice parameter *a*          real\*8 | real\*4   no          
  b          lattice parameter *b*          real\*8 | real\*4   no          
  c          lattice parameter *c*          real\*8 | real\*4   no          
  alpha      lattice angle *alpha*          real\*8 | real\*4   no          
  beta       lattice angle *beta*           real\*8 | real\*4   no          
  gamma      lattice angle *gamma*          real\*8 | real\*4   no          
  lenunits   units for lattice parameters   character(len=\*)   yes        angstrom
  angunits   units for lattice angles       character(len=\*)   yes        degree
  id         unique id                      character(len=\*)   yes         
  title      tag description                character(len=\*)   yes         
  dictref    dictionary reference           character(len=\*)   yes         
  fmt        format for reals               character(len=\*)   yes        f8.3

#### [cmlAddMetadata(xf, name, content, conv)]()

  argument   role                     type                optional   default
  ---------- ------------------------ ------------------- ---------- ---------
  xf         xml filehandle           type(xmlf\_t)       no          
  name       name (eg author)         character(len=\*)   no          
  content    value (eg Jon Wakelin)   character(len=\*)   no          
  conv       a convention             character(len=\*)   yes         

#### [cmlAddLattice(xf, cell, units, title, id, dictref, conv, lattType, spaceType)]()

  argument    role                             type                optional   default
  ----------- -------------------------------- ------------------- ---------- ---------
  xf          xml filehandle                   type(xmlf\_t)       no          
  cell(3,3)   the lattice vectors              real\*8 | real\*4   no          
  units       units for lattice vectors        character(len=\*)   yes         
  title       tag description                  character(len=\*)   yes         
  id          unique id                        character(len=\*)   yes         
  dictref     dictionary reference             character(len=\*)   yes         
  conv        a convention                     character(len=\*)   yes        f8.3
  lattType    lattice type (primitive |full)   character(len=\*)   yes         
  spaceType   space type (real | reciprocal)   character(len=\*)   yes         

#### [cmlAddProperty(xf, property, nvalues, ncols, nrows, dim, id, title, conv, dictref, ref, units, fmt)]()

  argument   role                                                                                        type                                             optional                 default
  ---------- ------------------------------------------------------------------------------------------- ------------------------------------------------ ------------------------ ---------
  xf         xml filehandle                                                                              type(xmlf\_t)                                    no                        
  property   a scalar, array or matrix property. Yes it takes all three, but see next three arguments.   charcter(len=\*) | integer | real\*8 | real\*4   no                        
  nvalue     number of value for vector                                                                  integer                                          no if writing a vector    
  ncols      number of rows for matrix                                                                   integer                                          no if writing a matrix    
  nrows      number of rows for matrix                                                                   integer                                          no if writing a matrix    
  dim        number of the fastest growing dimension in the matrix                                       integer                                          yes                       
  id         unique id                                                                                   character(len=\*)                                yes                       
  title      tag description                                                                             character(len=\*)                                yes                       
  conv       a convention                                                                                character(len=\*)                                yes                       
  dictref    dictionary reference                                                                        character(len=\*)                                yes                       
  ref        refernce to an atom (via the atom 'id' attribute)                                           character(len=\*)                                yes                       
  units      units for the property                                                                      character(len=\*)                                yes                       
  fmt        format for reals                                                                            character(len=\*)                                yes                      f8.3

#### [cmlAddParameter(xf, value, name, id, title, conv, cons, dictref, ref, role, units, fmt)]()

  argument   role                                                type                                                       optional   default
  ---------- --------------------------------------------------- ---------------------------------------------------------- ---------- ---------
  xf         xml filehandle                                      type(xmlf\_t)                                              no          
  value      the parameter                                       charcter(len=\*) | integer | real\*8 | real\*4 | logical   no          
  id         unique id                                           character(len=\*)                                          yes         
  title      tag description                                     character(len=\*)                                          yes         
  conv       a convention                                        character(len=\*)                                          yes         
  cons       a constraint                                        character(len=\*)                                          yes         
  role       the role of the parameter                           character(len=\*)                                          yes         
  dictref    dictionary reference                                character(len=\*)                                          yes         
  ref        refernce to an atom (via the atom 'id' attribute)   character(len=\*)                                          yes         
  units      units for the parameter                             character(len=\*)                                          yes         
  fmt        format for reals                                    character(len=\*)                                          yes        f8.3


