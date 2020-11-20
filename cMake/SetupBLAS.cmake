function(ConfigureBLAS TARGET)
  find_package(BLAS REQUIRED)
  if(NOT BLAS_FOUND)
    message(FATAL_ERROR
      "===============\n"
      "BLAS not found.\n"
      "===============\n")
  endif(NOT BLAS_FOUND)
  TARGET_INCLUDE_DIRECTORIES(${TARGET} PUBLIC ${BLAS_INCLUDE_PATH})
  TARGET_LINK_LIBRARIES(${TARGET} PUBLIC ${BLAS_LIBRARIES})
endfunction(ConfigureBLAS)
