include(FindPackageHandleStandardArgs)

find_library(METIS_LIBRARY NAMES metis)
find_path(METIS_INCLUDE_DIR NAMES metis.h)

find_package_handle_standard_args(METIS DEFAULT_MSG METIS_LIBRARY METIS_INCLUDE_DIR)

if(METIS_FOUND)
    mark_as_advanced(METIS_INCLUDE_DIR)
    mark_as_advanced(METIS_LIBRARY)
endif()

IF(METIS_FOUND AND NOT TARGET METIS::METIS)
    add_library(METIS::METIS IMPORTED SHARED)
    set_property(TARGET METIS::METIS PROPERTY IMPORTED_LOCATION ${METIS_LIBRARY})
    target_include_directories(METIS::METIS INTERFACE ${METIS_INCLUDE_DIR})
endif()
