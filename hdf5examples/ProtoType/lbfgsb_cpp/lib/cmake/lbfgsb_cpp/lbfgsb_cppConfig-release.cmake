#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "lbfgsb_cpp::lbfgsb_cpp" for configuration "Release"
set_property(TARGET lbfgsb_cpp::lbfgsb_cpp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(lbfgsb_cpp::lbfgsb_cpp PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/liblbfgsb_cpp.so"
  IMPORTED_SONAME_RELEASE "liblbfgsb_cpp.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS lbfgsb_cpp::lbfgsb_cpp )
list(APPEND _IMPORT_CHECK_FILES_FOR_lbfgsb_cpp::lbfgsb_cpp "${_IMPORT_PREFIX}/lib/liblbfgsb_cpp.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)