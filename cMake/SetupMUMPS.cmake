function(ConfigureMUMPS TARGET)

# Find mumps headers:
FIND_PATH(MUMPS_INCLUDE_PATH
    NAMES cmumps_c.h
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/include"
    "${PETSC_PATH}/arch-linux2-c-opt/include"
    "${PETSC_PATH}/arch-darwin-c-opt/include"
    "/usr/include"
    )

if(MUMPS_INCLUDE_PATH)
    message(STATUS "Mumps headers found at " ${MUMPS_INCLUDE_PATH})
else()
    message(STATUS "Mumps headers not found")
endif()

# Find mumps library:
FIND_LIBRARY(MUMPS_LIBRARIES
    NAMES cmumps
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/lib"
    "${PETSC_PATH}/arch-linux2-c-opt/lib"
    "${PETSC_PATH}/arch-darwin-c-opt/lib"
    "/usr/lib/x86_64-linux-gnu"
    "/usr/lib/aarch64-linux-gnu"
    )

if(MUMPS_LIBRARIES)
    message(STATUS "Mumps library found at " ${MUMPS_LIBRARIES})
else()
    message(STATUS "Mumps library not found")
endif()

if(MUMPS_INCLUDE_PATH AND MUMPS_LIBRARIES)
    SET(MUMPS_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${MUMPS_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${MUMPS_LIBRARIES})
endif()

endfunction(ConfigureMUMPS)

