// Tube thickness [m]:
thtube = 1e-3;
// Tube height [m]:
htube = 210e-3;
// Tube inner radius [m]:
rintube = 10e-3;
// Radius of the air domain around the tube [m]:
rinf = 200e-3;

// Number of mesh nodes on the tube:
nx = 20;
ny = 100;
// Characteristic mesh at infinity:
tc = 2.5e-2;
tcintube = 3e-3;


// Define the tube rectangle with a transfinite (i.e. structured) mesh:
Point(1) = {rintube, -htube/2, 0, 1};
Point(2) = {rintube+thtube, -htube/2, 0, 1};
Point(3) = {rintube+thtube, htube/2, 0, 1};
Point(4) = {rintube, htube/2, 0, 1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1}; 

Transfinite Line(1) = nx;
Transfinite Line(2) = ny;
Transfinite Line(3) = nx;
Transfinite Line(4) = ny;

Transfinite Surface(1);
Recombine Surface(1);

// Define the air around:
Point(5) = {0,-rinf,0, tc};
Point(6) = {0,rinf,0, tc};
Point(7) = {rinf,-rinf,0, tc};
Point(8) = {rinf,rinf,0, tc};
Point(9) = {0,-htube/2,0, tcintube};
Point(10) = {0,htube/2,0, tcintube};

Line(5) = {5,7};
Line(6) = {7,8};
Line(7) = {8,6};
Line(8) = {6,10};
Line(9) = {10,9};
Line(10) = {9,5};

Line Loop(2) = {6,7,8,9,10,5};
Plane Surface(2) = {1,2};

// Recombine triangles into quadrangles:
Recombine Surface {2};



// Define the physical regions:
tube = 1;
air = 2;
domainskin = 3;

Physical Surface(tube) = {1};
Physical Surface(air) = {2};
Physical Line(domainskin) = {5,6,7,8,9,10};


