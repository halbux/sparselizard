# Locate METIS include paths and libraries
# METIS can be found at http://glaros.dtc.umn.edu/gkhome/metis/metis/download

# This module defines
# METIS_INCLUDE_DIRS, where to find METIS.h, etc.
# METIS_LIBRARIES, the libraries to link against to use METIS.
# METIS_FOUND, If false, don't try to use METIS.

function(ConfigureMETIS TARGET)
  message(STATUS "Looking for METIS")

  # FIXME: windows paths were not tested
  FIND_PATH(METIS_INCLUDE_PATH
    NAMES metis.h
    PATHS
    ${CMAKE_INCLUDE_PATH}
    "[HKEY_CURRENT_USER\\software\\METIS\\include]"
    "[HKEY_CURRENT_USER\\METIS\\include]"
    $ENV{METISROOT}/include/
    $ENV{METIS_DIR}/include/
    /usr/local/include/
    /usr/include/
    )

  # FIXME: windows paths were not tested
  FIND_LIBRARY(METIS_LIBRARIES
    NAMES libmetis.so
    PATHS
    ${CMAKE_LIBRARY_PATH}
    "[HKEY_CURRENT_USER\\software\\METIS\\lib\\METIS]"
    "[HKEY_CURRENT_USER\\METIS\\lib\\METIS]"
    $ENV{METISROOT}/lib/METIS
    $ENV{METIS_DIR}/lib/
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    )

  IF(METIS_INCLUDE_PATH AND METIS_LIBRARIES)
    message(STATUS "Found METIS: ${METIS_LIBRARIES}")
    SET(METIS_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${METIS_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${METIS_LIBRARIES})
  ELSE(METIS_INCLUDE_PATH AND METIS_LIBRARIES)
    message(WARNING "METIS not found")
  ENDIF(METIS_INCLUDE_PATH AND METIS_LIBRARIES)


endfunction(ConfigureMETIS)
