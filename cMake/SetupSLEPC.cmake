# Locate SLEPC include paths and libraries
# SLEPC can be found at https://gitlab.com/slepc/slepc

# This module defines
# SLEPC_INCLUDE_DIRS, where to find SLEPC.h, etc.
# SLEPC_LIBRARIES, the libraries to link against to use SLEPC.
# SLEPC_FOUND, If false, don't try to use SLEPC.

function(ConfigureSLEPC TARGET)
  message(STATUS "Looking for SLEPC")

  # FIXME: windows paths were not tested
  FIND_PATH(SLEPC_INCLUDE_PATH
    NAMES slepc.h
    PATHS
    ${CMAKE_INCLUDE_PATH}
    "[HKEY_CURRENT_USER\\software\\slepc\\include]"
    "[HKEY_CURRENT_USER\\slepc\\include]"
    $ENV{SLEPCROOT}/include/
    $ENV{SLEPC_DIR}/include/
    /usr/local/include/
    /usr/include/
    /usr/include/slepc/
    REQUIRED ON)

  # FIXME: windows paths were not tested
  FIND_LIBRARY(SLEPC_LIBRARIES
    NAMES libslepc.so
    PATHS
    ${CMAKE_LIBRARY_PATH}
    "[HKEY_CURRENT_USER\\software\\slepc\\lib\\slepc]"
    "[HKEY_CURRENT_USER\\SLEPC\\lib\\slepc]"
    $ENV{SLEPCROOT}/lib/slepc
    $ENV{SLEPC_DIR}/lib/
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    REQUIRED ON)

  SET(SLEPC_FOUND YES)
  message(STATUS "Found SLEPC: ${SLEPC_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${SLEPC_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${SLEPC_LIBRARIES})

endfunction(ConfigureSLEPC)
