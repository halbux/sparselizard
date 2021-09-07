function(ConfigureMETIS TARGET)

# Find metis headers:
FIND_PATH(METIS_INCLUDE_PATH
    NAMES metis.h
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/include"
    "${PETSC_PATH}/arch-linux2-c-opt/include"
    "${PETSC_PATH}/arch-darwin-c-opt/include"
    )

if(METIS_INCLUDE_PATH)
    message(STATUS "Metis headers found at " ${METIS_INCLUDE_PATH})
else()
    message(STATUS "Metis headers not found")
endif()

# Find metis library:
FIND_LIBRARY(METIS_LIBRARIES
    NAMES metis
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/lib"
    "${PETSC_PATH}/arch-linux2-c-opt/lib"
    "${PETSC_PATH}/arch-darwin-c-opt/lib"
    )

if(METIS_LIBRARIES)
    message(STATUS "Metis library found at " ${METIS_LIBRARIES})
else()
    message(STATUS "Metis library not found")
endif()

if(METIS_INCLUDE_PATH AND METIS_LIBRARIES)
    SET(METIS_FOUND YES PARENT_SCOPE)

    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${METIS_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${METIS_LIBRARIES})
endif()

endfunction(ConfigureMETIS)