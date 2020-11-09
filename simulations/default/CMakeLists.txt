custom_add_executable_from_dir(default ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(default sparselizard)
custom_copy_file(default ${CMAKE_CURRENT_SOURCE_DIR} "*.msh;*.geo")
