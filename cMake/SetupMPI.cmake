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
    "/usr/include/x86_64-linux-gnu/mpi"
    "/usr/include/aarch64-linux-gnu/mpi"
    NO_DEFAULT_PATH
    )

if(MPI_INCLUDE_PATH)
    message(STATUS "MPI headers found at " ${MPI_INCLUDE_PATH})
else()
    message(STATUS "MPI headers not found (optional)")
endif()

# Find mpi library:
FIND_LIBRARY(MPI_LIBRARIES
    NAMES mpi
    PATHS
    "${MPI_PATH}/lib"
    "/usr/lib/x86_64-linux-gnu/openmpi/lib"
    "/usr/lib/aarch64-linux-gnu/openmpi/lib"
    NO_DEFAULT_PATH
    )

if(MPI_LIBRARIES)
    message(STATUS "MPI library found at " ${MPI_LIBRARIES})
else()
    message(STATUS "MPI library not found (optional)")
endif()

if(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
    SET(MPI_FOUND YES PARENT_SCOPE)
    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${MPI_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${MPI_LIBRARIES})
endif()

endfunction(ConfigureMPI)
