# Add library for dirs
function(custom_add_library_from_dir TARGET DIRLIST)
  foreach(d IN LISTS DIRLIST)
    file(GLOB SRC "${d}/*.cpp" "${d}/*.h" "${d}/*.hpp")
    list(APPEND TARGET_SRC ${SRC})
  endforeach()
  add_library(${TARGET} SHARED ${TARGET_SRC})
endfunction(custom_add_library_from_dir)

# Add executable for dirs
function(custom_add_executable_from_dir TARGET DIRLIST)
  foreach(d IN LISTS DIRLIST)
    file(GLOB SRC "${d}/*.cpp" "${d}/*.h" "${d}/*.hpp")
    list(APPEND TARGET_SRC ${SRC})
  endforeach()
  add_executable(${TARGET} ${TARGET_SRC})
endfunction(custom_add_executable_from_dir)

# Copy other files
function(custom_copy_file TARGET DIR GLOBS)
  foreach(g IN LISTS GLOBS)
    # message(FATAL_ERROR "${DIR}/${g}")
    file(GLOB SRC "${DIR}/${g}")
    foreach(f IN LISTS SRC)
      string(REGEX REPLACE "^.*/" "" DEST ${f})
      # TODO: maybe simlink instead of copy
      configure_file(${f} ${CMAKE_CURRENT_BINARY_DIR}/${DEST} COPYONLY)
    endforeach()
  endforeach()
endfunction(custom_copy_file)

# add tests
function(custom_add_test_from_dir TARGET LIBRARY DIRS)
  custom_add_executable_from_dir(${TARGET} ${DIRS})
  target_link_libraries(${TARGET} ${LIBRARY})
  add_test(${TARGET} ${TARGET})
endfunction()
