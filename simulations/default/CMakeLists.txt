custom_add_executable_from_dir(default ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(default sparselizard)
custom_symlink_file(default ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} "*.msh;*.geo;*.txt;*.csv;*.nas;*.slz")

