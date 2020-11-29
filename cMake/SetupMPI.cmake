function(ConfigureMPI TARGET)


# Find mpi headers:
FIND_PATH(MPI_INCLUDE_PATH
    NAMES mpi.h
    PATHS
    /usr/include/
    /usr/local/include/
    /usr/include/openmpi-x86_64
    /usr/include/mpich-x86_64
    )

if(MPI_INCLUDE_PATH)
    message(STATUS "MPI headers found at " ${MPI_INCLUDE_PATH})
else()
    message(STATUS "MPI HEADERS NOT FOUND")
endif()


# Find mpi library:
FIND_LIBRARY(MPI_LIBRARIES
    NAMES mpi
    PATHS
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    /usr/lib64
    /usr/lib64/openmpi/lib
    /usr/lib64/mpich/lib
    )

if(MPI_LIBRARIES)
    message(STATUS "MPI library found at " ${MPI_LIBRARIES})
else()
    message(STATUS "MPI LIBRARY NOT FOUND")
endif()


if(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
    SET(MPI_FOUND YES PARENT_SCOPE)

    TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${MPI_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${MPI_LIBRARIES})
endif()


endfunction(ConfigureMPI)

