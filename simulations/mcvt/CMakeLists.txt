custom_add_executable_from_dir(mcvt ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mcvt sparselizard)
custom_symlink_file(mcvt ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} "*.msh;*.geo;*.txt;*.csv;*.nas;*.slz")

