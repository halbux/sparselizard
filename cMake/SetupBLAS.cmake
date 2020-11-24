function(ConfigureBLAS TARGET)


# Find blas headers:
message(STATUS "Looking for blas headers at " ${PETSC_PATH})

FIND_PATH(BLAS_INCLUDE_PATH
    NAMES cblas.h
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/include"
    "${PETSC_PATH}/arch-linux2-c-opt/include"
    "${PETSC_PATH}/arch-darwin-c-opt/include"
    REQUIRED ON)

if(BLAS_INCLUDE_PATH)
    message(STATUS "Blas header cblas.h found at " ${BLAS_INCLUDE_PATH})
else()
    message(STATUS "BLAS HEADER CBLAS.H NOT FOUND")
endif()


# Find metis library:
message(STATUS "Looking for metis library at " ${PETSC_PATH})

FIND_LIBRARY(BLAS_LIBRARIES
    NAMES libopenblas.so
    PATHS
    "${PETSC_PATH}/arch-linux-c-opt/lib"
    "${PETSC_PATH}/arch-linux2-c-opt/lib"
    "${PETSC_PATH}/arch-darwin-c-opt/lib"
    REQUIRED ON)

if(BLAS_LIBRARIES)
    message(STATUS "Blas library found at " ${BLAS_LIBRARIES})
else()
    message(STATUS "BLAS LIBRARY NOT FOUND")
endif()


if(BLAS_INCLUDE_PATH AND BLAS_LIBRARIES)
    SET(BLAS_FOUND YES PARENT_SCOPE)
endif()


TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${BLAS_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${BLAS_LIBRARIES})


endfunction(ConfigureBLAS)

