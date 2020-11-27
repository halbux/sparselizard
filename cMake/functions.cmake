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
function(custom_copy_file TARGET FROMDIRS TODIR GLOBS)
    foreach(d IN LISTS FROMDIRS)
      foreach(g IN LISTS GLOBS)
        file(GLOB SRC "${d}/${g}")
        foreach(f IN LISTS SRC)
          string(REGEX REPLACE "^.*/" "" DEST ${f})
          configure_file(${f} ${TODIR}/${DEST} COPYONLY)
        endforeach()
      endforeach()
    endforeach()
endfunction(custom_copy_file)

