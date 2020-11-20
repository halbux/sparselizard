cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

project(Sparselizard LANGUAGES CXX VERSION 0.1)

set(DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Append to debug version
set(CMAKE_DEBUG_POSTFIX "_debug")

# Installation definitions
include(GNUInstallDirs)
# Aux functions
include(cMake/functions.cmake)
include(cMake/SetupBLAS.cmake)
include(cMake/SetupMPI.cmake)
include(cMake/SetupPETSC.cmake)
include(cMake/SetupSLEPC.cmake)
include(cMake/SetupGMSH.cmake)
include(cMake/SetupMETIS.cmake)

# add libsparselizard target
add_subdirectory(src)

# choose what to build
option(BUILD_EXAMPLES "Build examples instead of simulation." OFF)

if(BUILD_EXAMPLES)
  # add examples targets
  add_subdirectory(examples)
else(BUILD_EXAMPLES)
  # add simulations targets
  add_subdirectory(simulations)
endif(BUILD_EXAMPLES)
