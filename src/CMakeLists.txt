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

target_compile_features(sparselizard PRIVATE cxx_std_20)

ConfigureBLAS(sparselizard)
ConfigureMPI(sparselizard)
ConfigurePETSC(sparselizard)
ConfigureSLEPC(sparselizard)
ConfigureGMSH(sparselizard)
ConfigureMETIS(sparselizard)

# set configure file
if(${GMSH_FOUND})
    set(SPARSELIZARD_HAS_GMSH YES)
endif(${GMSH_FOUND})
if(${METIS_FOUND})
    set(SPARSELIZARD_HAS_METIS YES)
endif(${METIS_FOUND})
configure_file(config.h.in config.h)
target_include_directories(sparselizard PUBLIC ${CMAKE_CURRENT_BINARY_DIR})


# compiler warnings
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(sparselizard PRIVATE -Wall -Wextra -Wunreachable-code -Wpedantic)
endif()
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(sparselizard PRIVATE -Wweak-vtables -Wexit-time-destructors -Wglobal-constructors -Wmissing-noreturn )
endif()
if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(sparselizard PRIVATE /W4 /w44265 /w44061 /w44062 )
endif()
