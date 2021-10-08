# Thanks to S. Matsievskiy for bringing cmake to the project.

cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

project(Sparselizard LANGUAGES CXX)

set(DEFAULT_BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build" FORCE)

# Place library in build folder:
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -Og -Wall -Wextra -Wpedantic")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(NOT BLA_VENDOR)
    set(BLA_VENDOR CACHE STRING OpenBLAS)
endif()

if(NOT PETSC_DIR AND NOT(APPLE))
    set(PETSC_DIR CACHE STRING /usr/lib/petscdir)
endif()

# Add cmake packages
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_PREFIX_PATH $<$<BOOL:${APPLE}>:/opt/homebrew/opt/openblas> ${PETSC_DIR})

# Installation definitions
include(GNUInstallDirs)

# this is not the cmake-official way but it's a lot easier
set(INCLUDE_DIRS
    src
    src/mesh
    src/gausspoint
    src/io
    src/io/paraview
    src/io/nastran
    src/io/gmsh
    src/field
    src/resolution
    src/geometry
    src/expression
    src/expression/operation
    src/shapefunction
    src/shapefunction/hierarchical
    src/shapefunction/hierarchical/h1
    src/shapefunction/hierarchical/hcurl
    src/shapefunction/lagrange
    src/formulation
)
file(GLOB_RECURSE SRC_FILES "src/*.cpp")
file(GLOB_RECURSE SRC_HEADERS "src/*.h")
add_library(sparselizard SHARED ${SRC_FILES})
# Note: the following would allow hierarchical install, but sparselizard acts as if it's flattened
# target_sources(sparselizard PUBLIC FILE_SET HEADERS FILES ${SRC_HEADERS} BASE_DIRS src)
target_precompile_headers(sparselizard PUBLIC ${SRC_HEADERS})
target_include_directories(sparselizard PUBLIC ${INCLUDE_DIRS})
set_target_properties(sparselizard PROPERTIES PUBLIC_HEADER "${SRC_HEADERS}")

target_compile_options(sparselizard PUBLIC -fPIC -O3)

find_package(OpenMP REQUIRED)
target_link_libraries(sparselizard PUBLIC OpenMP::OpenMP_CXX)

find_package(MPI REQUIRED)
target_link_libraries(sparselizard PUBLIC MPI::MPI_CXX)

find_package(BLAS REQUIRED)
target_link_libraries(sparselizard PUBLIC BLAS::BLAS)

find_package(PETSC REQUIRED)
target_link_libraries(sparselizard PUBLIC PETSC::PETSC)

find_package(SLEPC)
target_link_libraries(sparselizard PUBLIC SLEPC::SLEPC)

find_package(GMSH)
if(GMSH_FOUND)
    target_link_libraries(sparselizard PUBLIC GMSH::GMSH)
    target_compile_definitions(sparselizard PRIVATE -DHAVE_GMSH)
endif()

find_package(METIS)

if(METIS_FOUND)
    target_link_libraries(sparselizard PUBLIC METIS::METIS)
    target_compile_definitions(sparselizard PRIVATE -DHAVE_METIS)
endif()

target_link_libraries(sparselizard PUBLIC cmumps)
target_compile_definitions(sparselizard PRIVATE
    -DHAVE_BLAS
    -DHAVE_MPI
    -DHAVE_MUMPS
    -DHAVE_PETSC
    -DHAVE_SLEPC
)

install(TARGETS sparselizard LIBRARY PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/sparselizard)

set(BUILD_SIMULATIONS YES)

if(BUILD_SIMULATIONS)
    add_subdirectory(simulations)
    message(STATUS "Building simulations")
endif()
