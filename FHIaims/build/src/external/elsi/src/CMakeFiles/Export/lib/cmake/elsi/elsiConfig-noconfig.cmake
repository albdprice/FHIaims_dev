#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "elsi::elpa" for configuration ""
set_property(TARGET elsi::elpa APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(elsi::elpa PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libelpa.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS elsi::elpa )
list(APPEND _IMPORT_CHECK_FILES_FOR_elsi::elpa "${_IMPORT_PREFIX}/lib/libelpa.a" )

# Import target "elsi::MatrixSwitch" for configuration ""
set_property(TARGET elsi::MatrixSwitch APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(elsi::MatrixSwitch PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libMatrixSwitch.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS elsi::MatrixSwitch )
list(APPEND _IMPORT_CHECK_FILES_FOR_elsi::MatrixSwitch "${_IMPORT_PREFIX}/lib/libMatrixSwitch.a" )

# Import target "elsi::OMM" for configuration ""
set_property(TARGET elsi::OMM APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(elsi::OMM PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libOMM.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS elsi::OMM )
list(APPEND _IMPORT_CHECK_FILES_FOR_elsi::OMM "${_IMPORT_PREFIX}/lib/libOMM.a" )

# Import target "elsi::NTPoly" for configuration ""
set_property(TARGET elsi::NTPoly APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(elsi::NTPoly PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libNTPoly.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS elsi::NTPoly )
list(APPEND _IMPORT_CHECK_FILES_FOR_elsi::NTPoly "${_IMPORT_PREFIX}/lib/libNTPoly.a" )

# Import target "elsi::fortjson" for configuration ""
set_property(TARGET elsi::fortjson APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(elsi::fortjson PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libfortjson.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS elsi::fortjson )
list(APPEND _IMPORT_CHECK_FILES_FOR_elsi::fortjson "${_IMPORT_PREFIX}/lib/libfortjson.a" )

# Import target "elsi::elsi" for configuration ""
set_property(TARGET elsi::elsi APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(elsi::elsi PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libelsi.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS elsi::elsi )
list(APPEND _IMPORT_CHECK_FILES_FOR_elsi::elsi "${_IMPORT_PREFIX}/lib/libelsi.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
