function(ConfigureBLAS TARGET)

# Find blas headers:
FIND_PATH(BLAS_INCLUDE_PATH
    NAMES cblas.h
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/include"
    "${PETSC_PATH}/arch-linux2-c-opt/include"
    "${PETSC_PATH}/arch-darwin-c-opt/include"
    "/usr/include/x86_64-linux-gnu"
    )

if(BLAS_INCLUDE_PATH)
    message(STATUS "Blas header cblas.h found at " ${BLAS_INCLUDE_PATH})
else()
    message(STATUS "Blas header cblas.h not found")
endif()

# Find blas library:
FIND_LIBRARY(BLAS_LIBRARIES
    NAMES openblas
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/lib"
    "${PETSC_PATH}/arch-linux2-c-opt/lib"
    "${PETSC_PATH}/arch-darwin-c-opt/lib"
    "/usr/lib/x86_64-linux-gnu"
    "/usr/lib/aarch64-linux-gnu"
    )

if(BLAS_LIBRARIES)
    message(STATUS "Blas library found at " ${BLAS_LIBRARIES})
else()
    message(STATUS "Blas library not found")
endif()

if(BLAS_INCLUDE_PATH AND BLAS_LIBRARIES)
    SET(BLAS_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${BLAS_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${BLAS_LIBRARIES})
endif()

endfunction(ConfigureBLAS)

