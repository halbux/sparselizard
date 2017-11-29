// Membrane length is 50 microns:
l = 50e-6;
// Electrode length is 10 microns:
lelectrode = 10e-6;
// Pillar length is  10 microns:
lp = 10e-6;
// Fluid half-circle radius:
rfluid = 6*l;
// Insulator thickness is 600 nanometers:
hinsulator = 6e-7;
// Vacuum gap height is 800 nanometers:
hvacuumgap = 8e-7;
// Total height is 2 microns:
h = 2e-6;

// 'lc' sets the mesh finesse:
lc = 1e-6;
// 'lcfluid' sets the mesh finesse in the fluid:
lcfluid = 30e-6;


Point(1) = {0, 0, 0, lc};
Point(2) = {l, 0, 0, lc};
Point(3) = {l, h, 0, lc};
Point(4) = {0, h, 0, lc};
Point(5) = {lp, hinsulator, 0, lc};
Point(6) = {l-lp, hinsulator, 0, lc};
Point(7) = {l-lp, hinsulator+hvacuumgap, 0, lc};
Point(8) = {lp, hinsulator+hvacuumgap, 0, lc};
Point(9) = {0, hinsulator, 0, lc};
Point(10) = {0, hinsulator+hvacuumgap, 0, lc};
Point(11) = {l, hinsulator, 0, lc};
Point(12) = {l, hinsulator+hvacuumgap, 0, lc};
Point(13) = {l/2-lelectrode/2, h, 0, lc};
Point(14) = {l/2+lelectrode/2, h, 0, lc};
Point(15) = {-rfluid+l/2, h, 0, lcfluid};
Point(16) = {rfluid+l/2, h, 0, lcfluid};
Point(17) = {l/2, h, 0, lcfluid};


Line(1) = {1, 2};
Line(2) = {2, 11};
Line(3) = {11, 12};
Line(4) = {12, 3};
Line(5) = {3, 14};
Line(6) = {14, 13};
Line(7) = {13, 4};
Line(8) = {4, 10};
Line(9) = {10, 9};
Line(10) = {9, 1};
Line(11) = {9, 5};
Line(12) = {5, 6};
Line(13) = {6, 11};
Line(14) = {12, 7};
Line(15) = {7, 8};
Line(16) = {8, 10};
Line(17) = {5, 8};
Line(18) = {6, 7};
Circle(19) = {16, 17, 15};
Line(20) = {15, 4};
Line(21) = {3, 16};

Line Loop(19) = {1, 2, -13, -12, -11, 10};
Plane Surface(1) = {19};
Line Loop(21) = {9, 11, 17, 16};
Plane Surface(2) = {21};
Line Loop(23) = {15, -17, 12, 18};
Plane Surface(3) = {23};
Line Loop(25) = {14, -18, 13, 3};
Plane Surface(4) = {25};
Line Loop(27) = {5, 6, 7, 8, -16, -15, -14, 4};
Plane Surface(5) = {27};
Line Loop(29) = {7, -20, -19, -21, 5, 6};
Plane Surface(6) = {29};

// Recombine as much as possible the triangles into quadrangles:
Recombine Surface(1);
Recombine Surface(2);
Recombine Surface(3);
Recombine Surface(4);
Recombine Surface(5);
Recombine Surface(6);


// The insulator:
Physical Surface(1) = {1};
// The membrane support pillars:
Physical Surface(2) = {2, 4};
// The vacuum gap:
Physical Surface(3) = {3};
// The membrane:
Physical Surface(4) = {5};
// The fluid:
Physical Surface(5) = {6};


// The ground and clamp line:
Physical Line(6) = {1};
// The membrane electrode:
Physical Line(7) = {6};
// The fluid boundary:
Physical Line(8) = {19,20,21};


