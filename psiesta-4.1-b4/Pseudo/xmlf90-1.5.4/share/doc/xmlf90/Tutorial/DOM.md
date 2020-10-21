FDOM
====

The FDOM is a DOM level 1.0 implementation written in F95. There are two
"gotchas" that the Fortran programmer should be aware of, and I can't
stress either strongly enough. Firstly the DOM, like many programming
languages, starts counting from 0, and not 1, for this reason all do
loops should run from 0 to length - 1. Secondly, we while we can not
return an object in Fortran *per se*, we can return a pointer to an
arbitrary structure containing mulitple members, including substructures
(which is very much like returning an object). Therefore, you must use
the pointer syntax:

    program main

     use xmlf90_dom
     
     type(fnode), pointer :: object
     type(fnode), pointer :: myNode
     
     object => getParentNode(myNode)
     

and not

    program main

     use xmlf90_dom
     
     type(fnode) :: object
     type(fnode) :: myNode
     
     object = getParentNode(myNode)
     

I can anticipate that these two issues (over running lists by counting
beyond the last item and trying to assign a pointer with =) will be the
cause of most programming errors using FDOM. But fairly soon you get
used to it, if you already program in other languages, or already use
pointers in Fortran there shouldn't be any problem at all.

Here is a list of the the methods implemented:

\

(Generic) Node Interface
------------------------

-   getNodeName(*node*)
-   getNodevalue)
-   getNodeType(*node*)
-   hasChildNodes(*node*)
-   hasAttributes(*node*)
-   getParentNode(*node*)
-   getFirstChild(*node*)
-   getLastChild(*node*)
-   getNextSibling(*node*)
-   getPreviousSibling(*node*)
-   getOwnerDocument(*node*)
-   getAttributes(*node*)
-   getChildNodes(*node*)
-   setNodeValue(*node, value*)
-   appendChild(*node, newChild*)
-   removeChild(*node, oldChild*)
-   replaceChild(*node, newChild, oldChild*)
-   cloneNode(*node, \[deep\]*)
-   isSameNode(*node, node2*)
-   insertBefore(*node, newChild, refChild*)

element Node Interface
----------------------

-   getTagName(*element*)
-   getElementsByTagName(*elment, tag*)
-   getAttribute(*element, name*)
-   getAttributeNode(*element, name*)
-   setAttribute(*element, name, value*)
-   setAttributeNode(*element, newAttribute*)
-   removeAttribute(*element, name*)

document Node Interface
-----------------------

-   createTextNode(*text*)
-   createAttribute(*name*)
-   createElement(*name*)
-   createComment(*data*)
-   getElementsByTagName(*document, tag*)

attribute Node Interface
------------------------

-   getName(*attr*)
-   getValue(*attr*)
-   setValue(*attr, value*)

nodeList Interface
------------------

-   item(*nodeList, i*)
-   getLength(*nodeList*)

namedNodeMap Interface
----------------------

-   item(*namedNodeMap, i*)
-   getLength(*namedNodeMap*)
-   getNamedItem(*namedNodeMap, name*)
-   setNamedItem(*namedNodeMap, name*)
-   removeNamedItem(*namedNodeMap, name*)

\
\
\
A partial list of the interfaces of the methods implemented follows. For
a full listing, please see the code in subdirectory `dom` of the main
distribution.\

-   -----------------------------------------------------------------------------------------------------
      method                                   arguments                  returns               DOM Level
      ---------------------------------------- -------------------------- --------------------- -----------
      getNodeName(node)                        type(fnode) :: node        string                1.0

      getNodeValue(node)                       type(fnode) :: node        string                1.0

      getNodeType(node)                        type(fnode) :: node        (integer code)        1.0

      getParentNode(node)                      type(fnode) :: node        type(fnode)           1.0

      getFirstChild(node)                      type(fnode) :: node        type(fnode)           1.0

      getLastChild(node)                       type(fnode) :: node        type(fnode)           1.0

      getPreviousSibling(node)                 type(fnode) :: node        type(fnode)           1.0

      getNextSibling(node)                     type(fnode) :: node        type(fnode)           1.0

      getOwnerDocument(node)                   type(fnode) :: node        type(fnode)           1.0

      getAttributes(node)                      type(fnode) :: node        type(fnamedNodeMap)   1.0

      getChildNodes(node)                      type(fnode) :: node        type(fnodeList)       1.0

      getOwnerDocument(node)                   type(fnode) :: node        type(fnode)           1.0

      appendChild(node, newChild)              type(fnode) :: node\       type(fnode)           1.0
                                               type(fnode) :: newChild                          

      removeChild(node, oldChild)              type(fnode) :: node\       type(fnode)           1.0
                                               type(fnode) :: oldChild                          

      replaceChild(node, newChild, oldChild)   type(fnode) :: node\       type(fnode)           1.0
                                               type(fnode) :: newChild\                         
                                               type(fnode) :: oldChild                          

      replaceChild(node, refChild, oldChild)   type(fnode) :: node\       type(fnode)           1.0
                                               type(fnode) :: refChild\                         
                                               type(fnode) :: oldChild                          

      hasChildren(node)                        type(fnode) :: node        logical               1.0

      hasAttributes                            type(fnode) :: node        logical               2.0

      isSameNode(node, node2)                  type(fnode) :: node\       logical               3.0
                                               type(fnode) :: node2                             
      -----------------------------------------------------------------------------------------------------

\
\

