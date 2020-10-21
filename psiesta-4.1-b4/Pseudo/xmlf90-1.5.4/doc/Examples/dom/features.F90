program features

  use xmlf90_dom

  implicit none

  type(fnode), pointer :: myDoc

  myDoc => parsefile(DOM_DATA_DIR//"test.xml" ) !! , verbose=.true.)
  print *, getNumberofAllocatedNodes()
  call dumpTree(myDoc)
  call xmlize(myDoc,"features.xml")
  call destroyNode(myDoc)
end program features
