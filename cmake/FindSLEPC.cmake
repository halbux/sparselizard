include(FindPackageHandleStandardArgs)

find_library(SLEPC_LIBRARY NAMES slepc slepc_real)
find_path(SLEPC_INCLUDE_DIR NAMES slepcsys.h HINTS
  /usr/lib/slepcdir
  /usr/include/slepc)

find_path(SLEPC_INCLUDE_DIR_PRIVATE NAMES slepcversion.h HINTS
  /usr/lib/slepcdir
  /usr/include/slepc
)

find_package_handle_standard_args(SLEPC DEFAULT_MSG SLEPC_LIBRARY SLEPC_INCLUDE_DIR)

if(SLEPC_FOUND)
  mark_as_advanced(SLEPC_INCLUDE_DIR)
  mark_as_advanced(SLEPC_INCLUDE_DIR_PRIVATE)
  mark_as_advanced(SLEPC_LIBRARY)
endif()

IF(SLEPC_FOUND AND NOT TARGET SLEPC::SLEPC)
  add_library(SLEPC::SLEPC IMPORTED SHARED)
  set_property(TARGET SLEPC::SLEPC PROPERTY IMPORTED_LOCATION ${SLEPC_LIBRARY})
  target_include_directories(SLEPC::SLEPC INTERFACE ${SLEPC_INCLUDE_DIR})
endif()
