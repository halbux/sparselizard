function(ConfigureGMSH TARGET)

# Find gmsh headers:
FIND_PATH(GMSH_INCLUDE_PATH
    NAMES gmsh.h
    PATHS
    "${GMSH_PATH}/include"
    )

if(GMSH_INCLUDE_PATH)
    message(STATUS "Gmsh headers found at " ${GMSH_INCLUDE_PATH})
else()
    message(STATUS "Gmsh headers not found (optional)")
endif()

# Find gmsh library:
FIND_LIBRARY(GMSH_LIBRARIES
    NAMES gmsh
    PATHS
    "${GMSH_PATH}/lib"
    )

if(GMSH_LIBRARIES)
    message(STATUS "Gmsh library found at " ${GMSH_LIBRARIES})
else()
    message(STATUS "Gmsh library not found (optional)")
endif()

if(GMSH_INCLUDE_PATH AND GMSH_LIBRARIES)
    SET(GMSH_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${GMSH_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${GMSH_LIBRARIES})
endif()

endfunction(ConfigureGMSH)

