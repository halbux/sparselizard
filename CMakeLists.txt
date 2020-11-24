# Thanks to S. Matsievskiy for bringing cmake to the project.

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(Sparselizard LANGUAGES CXX)

set(DEFAULT_BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)

set(PETSC_PATH "~/SLlibs/petsc" CACHE STRING "Provide the path to the petsc folder.")
set(GMSH_PATH "~/SLlibs/gmsh" CACHE STRING "Provide the path to the gmsh folder.")


# Installation definitions
include(GNUInstallDirs)
# Aux functions
include(cMake/functions.cmake)
include(cMake/SetupBLAS.cmake)
include(cMake/SetupGMSH.cmake)
include(cMake/SetupMETIS.cmake)
include(cMake/SetupPETSC.cmake)
include(cMake/SetupSLEPC.cmake)

# add libsparselizard target
add_subdirectory(src)

# add simulations targets
add_subdirectory(simulations)

