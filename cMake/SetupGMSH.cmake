# Locate GMSH include paths and libraries
# GMSH can be found at https://gitlab.onelab.info/gmsh/gmsh/

# This module defines
# GMSH_INCLUDE_DIRS, where to find GMSH.h, etc.
# GMSH_LIBRARIES, the libraries to link against to use GMSH.
# GMSH_FOUND, If false, don't try to use GMSH.

function(ConfigureGMSH TARGET)
  message(STATUS "Looking for GMSH")

  # FIXME: windows paths were not tested
  FIND_PATH(GMSH_INCLUDE_PATH
    NAMES gmsh.h
    PATHS
    ${CMAKE_INCLUDE_PATH}
    "[HKEY_CURRENT_USER\\software\\GMSH\\include]"
    "[HKEY_CURRENT_USER\\GMSH\\include]"
    $ENV{GMSHROOT}/include/
    /usr/local/include/
    /usr/include/
    )

  # FIXME: windows paths were not tested
  FIND_LIBRARY(GMSH_LIBRARIES
    NAMES libgmsh.so
    PATHS
    ${CMAKE_LIBRARY_PATH}
    "[HKEY_CURRENT_USER\\software\\GMSH\\lib\\GMSH]"
    "[HKEY_CURRENT_USER\\GMSH\\lib\\GMSH]"
    $ENV{GMSHROOT}/lib/GMSH
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    )

  IF(GMSH_INCLUDE_PATH AND GMSH_LIBRARIES)
    message(STATUS "Found GMSH: ${GMSH_LIBRARIES}")
    SET(GMSH_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${GMSH_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${GMSH_LIBRARIES})
  ELSE(GMSH_INCLUDE_PATH AND GMSH_LIBRARIES)
    message(WARNING "GMSH not found")
  ENDIF(GMSH_INCLUDE_PATH AND GMSH_LIBRARIES)

endfunction(ConfigureGMSH)
