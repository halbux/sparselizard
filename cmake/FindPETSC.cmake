include(FindPackageHandleStandardArgs)

find_library(PETSC_LIBRARY NAMES petsc petsc_real)
find_path(PETSC_INCLUDE_DIR NAMES petsc.h HINTS
  /usr/include/petsc
  /usr/lib/petscdir
  $<$<BOOL:${APPLE}>:/opt/homebrew/opt/petsc/include >
)

find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_LIBRARY PETSC_INCLUDE_DIR)

if(PETSC_FOUND)
  mark_as_advanced(PETSC_INCLUDE_DIR)
  mark_as_advanced(PETSC_LIBRARY)
endif()

IF(PETSC_FOUND AND NOT TARGET PETSC::PETSC)
  add_library(PETSC::PETSC IMPORTED SHARED)
  set_property(TARGET PETSC::PETSC PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARY})
  target_include_directories(PETSC::PETSC INTERFACE ${PETSC_INCLUDE_DIR})
endif()
