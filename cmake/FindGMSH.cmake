include(FindPackageHandleStandardArgs)

find_library(GMSH_LIBRARY NAMES gmsh)
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h)

find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH_LIBRARY GMSH_INCLUDE_DIR)

if(GMSH_FOUND)
    mark_as_advanced(GMSH_INCLUDE_DIR)
    mark_as_advanced(GMSH_LIBRARY)
endif()

IF(GMSH_FOUND AND NOT TARGET GMSH::GMSH)
    add_library(GMSH::GMSH IMPORTED SHARED)
    set_property(TARGET GMSH::GMSH PROPERTY IMPORTED_LOCATION ${GMSH_LIBRARY})
    target_include_directories(GMSH::GMSH INTERFACE ${GMSH_INCLUDE_DIR})
endif()
