
add_executable(default main.cpp)
target_link_libraries(default PUBLIC sparselizard)
install(TARGETS default RUNTIME)
set(DATA disk.geo disk.msh)
install(FILES ${DATA} DESTINATION ${CMAKE_INSTALL_BINDIR})
