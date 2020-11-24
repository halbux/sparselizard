// Choose the mesh size:
c = 0.1;

// Number of extrude layers:
n = 5;

Point(1) = {0, 0, 0, c};
Point(2) = {1, 0, 0, c};
Point(3) = {0, -1, 0, c};
Point(4) = {0, 1, 0, c};
Point(5) = {-1, 0, 0, c};

Circle(1) = {5, 1, 4};
Circle(2) = {4, 1, 2};
Circle(3) = {2, 1, 3};
Circle(4) = {3, 1, 5};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Recombine Surface{6};

Extrude{0,0,0.1}{Surface{6}; Layers{n}; Recombine;}

Physical Volume(1) = {1};
Physical Surface(2) = {15,19,23,27};
Physical Surface(3) = {28};
Physical Line(4) = {8,9,10,11};
