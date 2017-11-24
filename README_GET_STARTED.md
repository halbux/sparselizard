Welcome, new sparselizard user!

Sparselizard (Copyright (C) 2017-2018 Alexandre Halbach and Christophe Geuzaine, University of Liege, Belgium)
is an open source C++ finite element library meant to be user friendly and decently fast & parallelised.
It can handle a rather general set of problems in 1D, 2D and 3D such as mechanical, acoustic, thermal, electric and electromagnetic 
problems. Multiphysics problems, nonlinear problems or nonlinear multiphysics problems can be simulated as well.
The problems can be readily solved in time with a time-stepping resolution or with a user-friendly multiharmonic
resolution method. In the latter case the steady-state solution of a time-periodic problem can be obtained 
in a single step, for linear as well as for general nonlinear problems.
The library comes with hierarchical high order shape functions so that high order interpolations can be 
used with an interpolation order adapted to every unknown field and geometrical region.

For now sparselizard has been successfully tested on Linux and Mac (but not on Windows).
Working examples can be found in the 'examples' folder.

The widely-used open-source GMSH meshing software (www.gmsh.info) is recommended to mesh the geometry
and generate the .msh file required in the finite element simulation. The result files output by 
sparselizard are in .pos format supported by GMSH.

We hope you appreciate this library and wish you all the best with it!


##### FOLLOW THESE STEPS TO SET UP SPARSELIZARD :

1. Download sparselizard on your computer from https://gitlab.onelab.info/halbux/sparselizard.git 
	(e.g. run 'git clone https://gitlab.onelab.info/halbux/sparselizard.git' in a terminal).
	Open the folder in a terminal for the next steps.
	
2. Install the standard gfortran, gcc and g++ compilers required in the next step.
	On Ubuntu linux these can be installed easily with
	- sudo apt-get install gfortran
	- sudo apt-get install gcc
	- sudo apt-get install g++
	if not already there.

3. In the downloaded folder go to folder 'install_external_libs' and run all scripts in the numbered order.
	This installs in your home directory all third party libraries that are required for sparselizard
	(i.e. the widely-used petsc, slepc, openblas and fftw libraries). This may take some time and requires
	a working internet connection.
	If this step fails or if you wish to do it yourself follow the steps in the bash scripts,
	it should be quite easy to follow. If you do so do not forget to change the library path accordingly
	in the makefile and in 'run_sparselizard.sh'.

4. Now we are going to actually compile sparselizard. This should also be straightforward: just run 'make'
	in the terminal (or 'make -j4' if you have 4 computing cores).
	
5. You're done! Follow the examples in the 'examples' folder to discover sparselizard! But let's first 
	run a very simple test:
	
- Get the GMSH binary for your operating system and copy it to the sparselizard folder.
- Mesh the 'circle.geo' geometry by running './gmsh circle.geo -3' (3 because it is a 3D problem) or
		with './gmsh circle.geo' to mesh graphically. This creates a 'circle.msh' file which contains the mesh.
- Run './run_sparselizard.sh' in the terminal. This runs the code in 'main.cpp' that was compiled at the 
		previous step. When you edit 'main.cpp' you have to run again step 4 (but it will be much faster this time). 
- The previous step has created the 'u.pos' output file, which gives the exaggerated displacement of the top 
        surface in the thin cylinder geometry when the sides are clamped and a volume force is applied downwards. 
		Open it with './gmsh u.pos'. 
		You don't see anything or it looks weird? Don't worry, this is just because the simulation was performed 
		using very few hexahedra in the mesh but with an order 3 interpolation! 
		To visualise high order interpolations in GMSH do this: 
		
    - double click in the middle of the window then select 'All view options' at the bottom of the box that appeared
    - go to the 'General' tab
    - tick the 'Adapt visualization grid' box
    - Set 'Maximum recursion level' to 3 and 'Target visualization error' to the smallest possible value then 
        press enter. Now you have a finer solution! Since the solution is a mechanical displacement you might 
        want to see the deflection in 3D by double clicking in the middle of the window then selecting 
        'View vector display' >> 'Displacement' with factor 1.
		
		
		
