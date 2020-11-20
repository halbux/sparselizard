# Locate PETSc include paths and libraries
# PETSc can be found at https://gitlab.com/petsc/petsc

# This module defines
# PETSC_INCLUDE_DIRS, where to find petsc.h, etc.
# PETSC_LIBRARIES, the libraries to link against to use petsc.
# PETSC_FOUND, If false, don't try to use petsc.

function(ConfigurePETSC TARGET)
  message(STATUS "Looking for PETSc")

  # FIXME: windows paths were not tested
  FIND_PATH(PETSC_INCLUDE_PATH
    NAMES petsc.h
    PATHS
    ${CMAKE_INCLUDE_PATH}
    "[HKEY_CURRENT_USER\\software\\petsc\\include]"
    "[HKEY_CURRENT_USER\\petsc\\include]"
    $ENV{PETSCROOT}/include/
    $ENV{PETSC_DIR}/include/
    /usr/local/include/
    /usr/include/
    /usr/include/petsc/
    REQUIRED ON)

  # FIXME: windows paths were not tested
  FIND_LIBRARY(PETSC_LIBRARIES
    NAMES libpetsc.so
    PATHS
    ${CMAKE_LIBRARY_PATH}
    "[HKEY_CURRENT_USER\\software\\petsc\\lib\\petsc]"
    "[HKEY_CURRENT_USER\\petsc\\lib\\petsc]"
    $ENV{PETSCROOT}/lib/petsc
    $ENV{PETSC_DIR}/lib/
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    REQUIRED ON)

  SET(PETSC_FOUND YES)
  message(STATUS "Found PETSc: ${PETSC_LIBRARIES}")
  TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${PETSC_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${PETSC_LIBRARIES})

endfunction(ConfigurePETSC)
