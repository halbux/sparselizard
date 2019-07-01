// Choose the mesh size:
c = 0.004;

// Truss height, length, thickness:
h = 1; l = 0.5; th = 0.01;

Point(1) = {0, 0, 0, c};
Point(2) = {th, 0, 0, c};
Point(3) = {0, h, 0, c};
Point(4) = {th, h+th, 0, c};
Point(5) = {-l, h, 0, c};
Point(6) = {-l, h+th, 0, c};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 6};
Line(4) = {6, 5};
Line(5) = {5, 3};
Line(6) = {3, 1};

Curve Loop(1) = {5, 6, 1, 2, 3, 4};
Plane Surface(1) = {1};

solid = 1; clamp = 2; load = 3;

Physical Surface(solid) = {1};
Physical Line(clamp) = {1};
Physical Line(load) = {4};
