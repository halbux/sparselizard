// Choose the mesh size:
c = 0.5e-4;

// Number of extrude layers:
n = 2;

// Radius:
R = 800e-6 / 2;
// Thickness:
th = 10e-6;

Point(1) = {0, 0, 0, c};
Point(2) = {R, 0, 0, c};
Point(3) = {0, -R, 0, c};
Point(4) = {0, R, 0, c};
Point(5) = {-R, 0, 0, c};

Circle(1) = {5, 1, 4};
Circle(2) = {4, 1, 2};
Circle(3) = {2, 1, 3};
Circle(4) = {3, 1, 5};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Recombine Surface{6};

Extrude{0,0,th/2}{Surface{6}; Layers{n}; Recombine;}
Extrude{0,0,th/2}{Surface{28}; Layers{n}; Recombine;}

Physical Volume(1) = {1};
Physical Volume(2) = {2};
Physical Surface(3) = {15,19,23,27};
Physical Surface(4) = {50};
