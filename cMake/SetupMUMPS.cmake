# Locate MUMPS include paths and libraries
# MUMPS can be found at https://salsa.debian.org/science-team/mumps

# This module defines
# MUMPS_INCLUDE_DIRS, where to find mumps.h, etc.
# MUMPS_LIBRARIES, the libraries to link against to use MUMPS.
# MUMPS_FOUND, If false, don't try to use MUMPS.

function(ConfigureMUMPS TARGET)
  message(STATUS "Looking for MUMPS")

  # FIXME: windows paths were not tested
  FIND_PATH(MUMPS_INCLUDE_PATH
    NAMES cmumps_c.h
    PATHS
    ${CMAKE_INCLUDE_PATH}
    "[HKEY_CURRENT_USER\\software\\MUMPS\\include]"
    "[HKEY_CURRENT_USER\\MUMPS\\include]"
    $ENV{MUMPSROOT}/include/
    $ENV{MUMPS_DIR}/include/
    /usr/local/include/
    /usr/include/
    )

  # FIXME: windows paths were not tested
  FIND_LIBRARY(MUMPS_LIBRARIES
    NAMES libmumps_common.so
    PATHS
    ${CMAKE_LIBRARY_PATH}
    "[HKEY_CURRENT_USER\\software\\MUMPS\\lib\\MUMPS]"
    "[HKEY_CURRENT_USER\\MUMPS\\lib\\MUMPS]"
    $ENV{MUMPSROOT}/lib/MUMPS
    $ENV{MUMPS_DIR}/lib/
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    )

  IF(MUMPS_INCLUDE_PATH AND MUMPS_LIBRARIES)
    message(STATUS "Found MUMPS: ${MUMPS_LIBRARIES}")
    SET(MUMPS_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${MUMPS_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${MUMPS_LIBRARIES})
  ELSE(MUMPS_INCLUDE_PATH AND MUMPS_LIBRARIES)
    message(FATAL_ERROR
      "================\n"
      "MUMPS not found.\n"
      "================\n")
  ENDIF(MUMPS_INCLUDE_PATH AND MUMPS_LIBRARIES)


endfunction(ConfigureMUMPS)
