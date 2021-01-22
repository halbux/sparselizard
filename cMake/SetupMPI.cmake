# Run the default example with MPICH using:
#
# module avail mpi
# module load mpi/mpich-x86_64
# mpirun -n 4 ./default

function(ConfigureMPI TARGET)


# Find mpi headers:
FIND_PATH(MPI_INCLUDE_PATH
    NAMES mpi.h
    PATHS
    "${MPI_PATH}/include"
    NO_DEFAULT_PATH
    )

if(MPI_INCLUDE_PATH)
    message(STATUS "MPI headers found at " ${MPI_INCLUDE_PATH})
else()
    message(STATUS "MPI HEADERS NOT FOUND (OPTIONAL)")
endif()


# Find mpi library:
FIND_LIBRARY(MPI_LIBRARIES
    NAMES mpi
    PATHS
    "${MPI_PATH}/lib"
    NO_DEFAULT_PATH
    )

if(MPI_LIBRARIES)
    message(STATUS "MPI library found at " ${MPI_LIBRARIES})
else()
    message(STATUS "MPI LIBRARY NOT FOUND (OPTIONAL)")
endif()


if(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
    SET(MPI_FOUND YES PARENT_SCOPE)

    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${MPI_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${MPI_LIBRARIES})
endif()


endfunction(ConfigureMPI)

