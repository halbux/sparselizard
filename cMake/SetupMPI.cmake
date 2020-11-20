function(ConfigureMPI TARGET)
  find_package(MPI REQUIRED)
  if(NOT MPI_FOUND)
    message(FATAL_ERROR
      "==============\n"
      "MPI not found.\n"
      "==============\n")
  endif(NOT MPI_FOUND)
  TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${MPI_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${MPI_LIBRARIES})
endfunction(ConfigureMPI)
