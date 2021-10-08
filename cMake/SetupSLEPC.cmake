function(ConfigureSLEPC TARGET)

# Find slepc headers:
FIND_PATH(SLEPC_INCLUDE_PATH
    NAMES slepc.h
    HINTS
    "${PETSC_PATH}/arch-linux-c-opt/include"
    "${PETSC_PATH}/arch-linux2-c-opt/include"
    "${PETSC_PATH}/arch-darwin-c-opt/include"
    "/usr/lib/slepc/include"
    NO_DEFAULT_PATH
    )

if(SLEPC_INCLUDE_PATH)
    message(STATUS "Slepc headers found at " ${SLEPC_INCLUDE_PATH})
else()
    message(STATUS "Slepc headers not found")
endif()

# Find slepc library:
FIND_LIBRARY(SLEPC_LIBRARIES
    NAMES slepc slepc_real
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/lib"
    "${PETSC_PATH}/arch-linux2-c-opt/lib"
    "${PETSC_PATH}/arch-darwin-c-opt/lib"
    "/usr/lib/slepc/lib"
    NO_DEFAULT_PATH
    )

FIND_LIBRARY(SLEPC_LIBRARIES NAMES slepc_real)

if(SLEPC_LIBRARIES)
    message(STATUS "Slepc library found at " ${SLEPC_LIBRARIES})
else()
    message(STATUS "Slepc library not found")
endif()

if(SLEPC_INCLUDE_PATH AND SLEPC_LIBRARIES)
    SET(SLEPC_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${SLEPC_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${SLEPC_LIBRARIES})
endif()

endfunction(ConfigureSLEPC)

