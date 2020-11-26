# Thanks to S. Matsievskiy for bringing cmake to the project.

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(Sparselizard LANGUAGES CXX)

set(DEFAULT_BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)

set(PETSC_PATH "~/SLlibs/petsc" CACHE STRING "Provide the path to the petsc folder.")
set(GMSH_PATH "~/SLlibs/gmsh" CACHE STRING "Provide the path to the gmsh folder.")

# Default flags:
SET(BLAS_FOUND NO)
SET(GMSH_FOUND NO)
SET(METIS_FOUND NO)
SET(MUMPS_FOUND NO)
SET(PETSC_FOUND NO)
SET(SLEPC_FOUND NO)

# Installation definitions
include(GNUInstallDirs)
# Aux functions
include(cMake/functions.cmake)
include(cMake/SetupBLAS.cmake)
include(cMake/SetupGMSH.cmake)
include(cMake/SetupMETIS.cmake)
include(cMake/SetupMUMPS.cmake)
include(cMake/SetupPETSC.cmake)
include(cMake/SetupSLEPC.cmake)

# Add libsparselizard target
add_subdirectory(src)

# Add simulations targets
add_subdirectory(simulations)

