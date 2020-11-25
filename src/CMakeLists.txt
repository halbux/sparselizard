set(SRC_DIR_LIST
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/resolution
  ${CMAKE_CURRENT_SOURCE_DIR}/shapefunction
  ${CMAKE_CURRENT_SOURCE_DIR}/shapefunction/lagrange
  ${CMAKE_CURRENT_SOURCE_DIR}/shapefunction/hierarchical
  ${CMAKE_CURRENT_SOURCE_DIR}/shapefunction/hierarchical/h1
  ${CMAKE_CURRENT_SOURCE_DIR}/shapefunction/hierarchical/hcurl
  ${CMAKE_CURRENT_SOURCE_DIR}/shapefunction/hierarchical/meca
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh
  ${CMAKE_CURRENT_SOURCE_DIR}/formulation
  ${CMAKE_CURRENT_SOURCE_DIR}/field
  ${CMAKE_CURRENT_SOURCE_DIR}/io
  ${CMAKE_CURRENT_SOURCE_DIR}/io/nastran
  ${CMAKE_CURRENT_SOURCE_DIR}/io/paraview
  ${CMAKE_CURRENT_SOURCE_DIR}/io/gmsh
  ${CMAKE_CURRENT_SOURCE_DIR}/geometry
  ${CMAKE_CURRENT_SOURCE_DIR}/expression
  ${CMAKE_CURRENT_SOURCE_DIR}/expression/operation
  ${CMAKE_CURRENT_SOURCE_DIR}/gausspoint)

custom_add_library_from_dir(sparselizard "${SRC_DIR_LIST}")

target_include_directories(sparselizard PUBLIC ${SRC_DIR_LIST})

target_compile_features(sparselizard PRIVATE cxx_std_11)
target_compile_options(sparselizard PUBLIC -fopenmp -fPIC -no-pie -O3)

ConfigureBLAS(sparselizard)
ConfigurePETSC(sparselizard)
ConfigureSLEPC(sparselizard)
ConfigureGMSH(sparselizard)
ConfigureMETIS(sparselizard)

if(NOT (${BLAS_FOUND} AND ${METIS_FOUND} AND ${PETSC_FOUND} AND ${SLEPC_FOUND}))
    message(FATAL_ERROR "At least one required library was not found.")
endif()

# set configure file
if(${BLAS_FOUND})
    set(SPARSELIZARD_HAS_BLAS YES)
endif()
if(${GMSH_FOUND})
    set(SPARSELIZARD_HAS_GMSH YES)
endif()
if(${METIS_FOUND})
    set(SPARSELIZARD_HAS_METIS YES)
endif()
if(${PETSC_FOUND})
    set(SPARSELIZARD_HAS_PETSC YES)
endif()
if(${SLEPC_FOUND})
    set(SPARSELIZARD_HAS_SLEPC YES)
endif()
configure_file(config.h.in config.h)
target_include_directories(sparselizard PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

