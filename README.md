# Build instructions

Regular build:
```bash
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

---

For building with debugging symbols use:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

For extra build output use:
```bash
cmake --build . -j$(nproc) -v
```

It may be convenient to use tool GUI:
```bash
cmake-gui -S . -B build
```
In *Advanced Mode* you may set paths for auxiliary libraries (e.g. BLAS).

Build all examples:
```bash
cmake .. -DBUILD_EXAMPLES=ON
```

Build specific example:
```bash
cmake .. -DBUILD_EXAMPLES=ON -DBUILD_ALL_EXAMPLES=OFF \
    -DBUILD_ADVECTION_DIFFUSION_STABILIZATION_SANDBOX_2D=ON
```

Each example folder has `CMakeLists.txt` file, where you can lookup build flag for the example.

# Add project

Simulation projects are located under `simulations`.
In order to create a new simulation:

1. Copy `simulations/default` folder with different name. Let's say that the new folder is
   `simulations/newsim`
1. Replace target name `default` with the unique one in `simulations/newsim/CMakeLists.txt`
1. Add line `add_subdirectory(newsim)` to the `simulations/CMakeLists.txt`
1. Build as usual. Executable file will be located in `build/simulations/newsim` folder
