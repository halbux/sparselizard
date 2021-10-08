# Build instructions

Sparselizard depends on mandatory and optional external libraries.

On debian Linux, run this command to install these libraries globally
```sudo apt-get install -y libopenblas-dev libmetis-dev libopenmpi-dev libmumps-dev petsc-dev slepc-dev libgmsh-dev```
or
run the scripts in the 'install_external_libs' folder to download and install these libraries locally.

Then:
```bash
mkdir build; cmake -B build .; cd build && make -j $(getconf _NPROCESSORS_ONLN)
```

---

Provide a custom path to the petsc, gmsh (optional) or mpi (optional) folder with:
```bash
cmake -B build . -DPETSC_PATH=/yourpath/petsc -DGMSH_PATH=/yourpath/gmsh -DMPI_PATH=/yourpath/mpi
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
1. Configure and build. Executable file will be located in `build/simulations/newsim` folder

