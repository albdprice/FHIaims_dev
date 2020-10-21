#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "aims::aims" for configuration ""
set_property(TARGET aims::aims APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(aims::aims PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/aims.190903.scalapack.mpi.x"
  )

list(APPEND _IMPORT_CHECK_TARGETS aims::aims )
list(APPEND _IMPORT_CHECK_FILES_FOR_aims::aims "${_IMPORT_PREFIX}/bin/aims.190903.scalapack.mpi.x" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
