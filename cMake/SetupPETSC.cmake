function(ConfigurePETSC TARGET)


# Find petsc headers:
FIND_PATH(PETSC_INCLUDE_PATH
    NAMES petsc.h
    PATHS
    "${PETSC_PATH}/include"
    NO_DEFAULT_PATH
    )

if(PETSC_INCLUDE_PATH)
    message(STATUS "Petsc header petsc.h found at " ${PETSC_INCLUDE_PATH})
else()
    message(STATUS "PETSC HEADER PETSC.H NOT FOUND")
endif()

FIND_PATH(PETSCCONF_INCLUDE_PATH
    NAMES petscconf.h
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/include"
    "${PETSC_PATH}/arch-linux2-c-opt/include"
    "${PETSC_PATH}/arch-darwin-c-opt/include"
    NO_DEFAULT_PATH
    )

if(PETSCCONF_INCLUDE_PATH)
    message(STATUS "Petsc header petscconf.h found at " ${PETSCCONF_INCLUDE_PATH})
else()
    message(STATUS "PETSC HEADER PETSCCONF.H NOT FOUND")
endif()


# Find petsc library:
FIND_LIBRARY(PETSC_LIBRARIES
    NAMES petsc
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/lib"
    "${PETSC_PATH}/arch-linux2-c-opt/lib"
    "${PETSC_PATH}/arch-darwin-c-opt/lib"
    NO_DEFAULT_PATH
    )

if(PETSC_LIBRARIES)
    message(STATUS "Petsc library found at " ${PETSC_LIBRARIES})
else()
    message(STATUS "PETSC LIBRARY NOT FOUND")
endif()


if(PETSC_INCLUDE_PATH AND PETSCCONF_INCLUDE_PATH AND PETSC_LIBRARIES)
    SET(PETSC_FOUND YES PARENT_SCOPE)

    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${PETSC_INCLUDE_PATH} ${PETSCCONF_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${PETSC_LIBRARIES})
endif()


endfunction(ConfigurePETSC)

