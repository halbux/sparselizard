function(ConfigureCATCH TARGET)
  find_package(Catch2 REQUIRED)
  if(NOT Catch2_FOUND)
    message(FATAL_ERROR
      "================\n"
      "CATCH not found.\n"
      "================\n")
  endif(NOT Catch2_FOUND)
  TARGET_LINK_LIBRARIES(${TARGET} Catch2::Catch2)
endfunction(ConfigureCATCH)
