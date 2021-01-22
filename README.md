# Build instructions

Run the scripts in the 'install_external_libs' folder then:
```bash
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

---

To provide a custom path to the petsc, gmsh (optional) or mpi (optional) folder use:
```bash
cmake .. -DPETSC_PATH=/yourpath/petsc -DGMSH_PATH=/yourpath/gmsh -DMPI_PATH=/yourpath/mpi
```

It may be convenient to use the cmake GUI:
```bash
cmake-gui
```

# Add project

Simulation projects are located under `simulations`.
In order to create a new simulation:

1. Copy `simulations/default` folder with different name. Let's say that the new folder is
   `simulations/newsim`
1. Replace target name `default` with the new one in `simulations/newsim/CMakeLists.txt`
1. Add line `add_subdirectory(newsim)` to `simulations/CMakeLists.txt`
1. Build as usual. Executable file will be located in `build/simulations/newsim` folder

