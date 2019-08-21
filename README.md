User-friendly, multiphysics, multithreaded, robust FEM c++ library routinely used for next-generation microtransducer design.

You will find all information you need at www.sparselizard.org (visit the forum if you have any question).

## Installing Dependencies
Sparselizard requires some external dependencies before compilation. Compiling the dependencies typically only needs to be done once, but may need updating on occasion. Below are a couple of different ways to get the dependencies.

### Linking a prior Sparselizard libs directory
If you already have an `~/SLlibs` set of external libs (or for some reason have multiple folders of the `sparselizard` source code) and would not like to redownload/compile, you can run the below command to create a softlink to the appropriate location. The folder in which you execute the `link.sh` script does not matter.

``` bash
./external_libs/scripts/link.sh path/to/SLlibs
```

You will now notice there is a `libs` softlink to `path/to/SLlibs`. You can remove this link by executing:

``` bash
unlink ./external_libs/libs
```

### Downloading a fresh copy of external libraries
You can download/compile the external libraries with the below command. The directory in which you invoke `install.sh` does not matter.

``` bash
./external_libs/scripts/install.sh
```
