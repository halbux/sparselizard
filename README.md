## Release status

[![latest released version(s)](https://repology.org/badge/latest-versions/sparselizard.svg)](https://repology.org/project/sparselizard/versions)
[![Release status](https://repology.org/badge/tiny-repos/sparselizard.svg)](https://repology.org/metapackage/sparselizard/versions)


# Installation 

## Linux 

### Arch Linux

Sparselizard is available in the AUR: (https://aur.archlinux.org/packages/sparselizard/).

### Ubuntu

Sparselizard is packaged in a PPA.
```bash
sudo add-apt-repository ppa:js-reynaud/sparselizard-v0
sudo apt-get update
sudo apt-get install sparslizard
```

### Build instructions

Run the scripts in the 'install_external_libs' folder then configure and build:
```bash
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

---

Provide a custom path to the petsc, gmsh (optional) or mpi (optional) folder with:
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
1. Configure and build. Executable file will be located in `build/simulations/newsim` folder

