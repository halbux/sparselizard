// Choose the mesh size:
c = 0.1;

Point(1) = {-1, -1, 0, c};
Point(2) = {1, -1, 0, 3*c};
Point(3) = {1, 1, 0, c};
Point(4) = {-1, 1, 0, 2*c};
Point(5) = {0, 0, 0, c};

Circle(1) = {2, 5, 3};
Circle(2) = {3, 5, 4};
Circle(3) = {4, 5, 1};
Circle(4) = {1, 5, 2};

Line(5) = {1, 5};
Line(6) = {5, 4};
Line(7) = {5, 3};
Line(8) = {5, 2};

Line Loop(1) = {3, 5, 6};
Plane Surface(1) = {1};
Line Loop(2) = {6, -2, -7};
Plane Surface(2) = {2};
Line Loop(3) = {1, -7, 8};
Plane Surface(3) = {3};
Line Loop(4) = {5, 8, -4};
Plane Surface(4) = {4};

Recombine Surface{1};
Recombine Surface{3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Line(5) = {1};
Physical Line(6) = {2};
Physical Line(7) = {3};
Physical Line(8) = {4};
