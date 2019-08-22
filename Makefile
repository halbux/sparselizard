mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(dir $(mkfile_path))
external_libs_dir := $(current_dir)external_libs/libs

##### THESE ARE THE REQUIRED LIBRARIES:

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
LIBS = \
	-L $(external_libs_dir)/petsc/arch-linux2-c-opt/lib \
	-l openblas \
	-l petsc \
	-l slepc
INCL = \
	-I $(external_libs_dir)/petsc/include/petsc/mpiuni \
	-I $(external_libs_dir)/petsc/arch-linux2-c-opt/externalpackages/git.openblas \
	-I $(external_libs_dir)/petsc/include/ \
	-I $(external_libs_dir)/petsc/arch-linux2-c-opt/include/
endif
ifeq ($(UNAME), Darwin)
LIBS = \
	-L $(external_libs_dir)/petsc/arch-darwin-c-opt/lib \
	-l openblas \
	-l petsc \
	-l slepc
INCL = \
	-I $(external_libs_dir)/petsc/include/petsc/mpiuni \
	-I $(external_libs_dir)/petsc/arch-darwin-c-opt/externalpackages/git.openblas \
	-I $(external_libs_dir)/petsc/include/ \
	-I $(external_libs_dir)/petsc/arch-darwin-c-opt/include/
endif


# $@ is the filename representing the target.
# $< is the filename of the first prerequisite.
# $^ the filenames of all the prerequisites.
# $(@D) is the file path of the target file. 
# D can be added to all of the above.

CXX = g++ -fopenmp # openmp is used for parallel sorting on Linux
CXX_FLAGS= -std=c++11 -O3 -Wno-return-type # Do not warn when not returning anything (e.g. in virtual functions)

# List of all directories containing the headers:
INCLUDES = -I src -I src/field -I src/expression -I src/expression/operation -I src/shapefunction -I src/formulation -I src/shapefunction/hierarchical -I src/shapefunction/hierarchical/h1 -I src/shapefunction/hierarchical/hcurl -I src/shapefunction/hierarchical/meca -I src/gausspoint -I src/shapefunction/lagrange -I src/mesh -I src/io -I src/io/gmsh -I src/io/paraview -I src/io/nastran -I src/resolution -I src/geometry
# List of all .cpp source files:
CPPS= $(wildcard src/*.cpp) $(wildcard src/field/*.cpp) $(wildcard src/expression/*.cpp) $(wildcard src/expression/operation/*.cpp) $(wildcard src/shapefunction/*.cpp) $(wildcard src/formulation/*.cpp) $(wildcard src/shapefunction/hierarchical/*.cpp) $(wildcard src/shapefunction/hierarchical/h1/*.cpp) $(wildcard src/shapefunction/hierarchical/meca/*.cpp) $(wildcard src/shapefunction/hierarchical/hcurl/*.cpp) $(wildcard src/gausspoint/*.cpp) $(wildcard src/shapefunction/lagrange/*.cpp) $(wildcard src/mesh/*.cpp) $(wildcard src/io/*.cpp) $(wildcard src/io/gmsh/*.cpp) $(wildcard src/io/paraview/*.cpp) $(wildcard src/io/nastran/*.cpp) $(wildcard src/resolution/*.cpp) $(wildcard src/geometry/*.cpp)
# Final binary name:
BIN = sparselizard
# Put all generated stuff to this build directory:
BUILD_DIR = ./build


# Same list as CPP but with the .o object extension:
OBJECTS=$(CPPS:%.cpp=$(BUILD_DIR)/%.o)
# Gcc/Clang will create these .d files containing dependencies.
DEP = $(OBJECTS:%.o=%.d)

all: $(OBJECTS)
	# The main is always recompiled (it could have been replaced):
	$(CXX) $(CXX_FLAGS) $(LIBS) $(INCL) $(INCLUDES) -c main.cpp -o $(BUILD_DIR)/main.o
	# Linking objects:
	$(CXX) $(BUILD_DIR)/main.o $(OBJECTS) $(LIBS) -o $(BIN)
	
# Include all .d files
-include $(DEP)

$(BUILD_DIR)/%.o: %.cpp
	# Create the folder of the current target in the build directory:
	mkdir -p $(@D)
	# Compile .cpp file. MMD creates the dependencies.
	$(CXX) $(CXX_FLAGS) $(LIBS) $(INCL) $(INCLUDES) -MMD -c $< -o $@
	

clean :
    # Removes all files created.
	rm -rf $(BUILD_DIR)
	rm -f $(BIN)
