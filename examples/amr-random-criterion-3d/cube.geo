// Choose the mesh size:
c = 0.2;

Point(1) = {0, 0, 0, c};
Point(2) = {1, 0, 0, c};
Point(3) = {1, 1, 0, c};
Point(4) = {0, 1, 0, c};
Point(5) = {0, 0, 1, c};
Point(6) = {1, 0, 1, c};
Point(7) = {1, 1, 1, c};
Point(8) = {0, 1, 1, c};

Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};
Line(5) = {1, 4};
Line(6) = {4, 8};
Line(7) = {8, 5};
Line(8) = {6, 7};
Line(9) = {7, 3};
Line(10) = {3, 2};
Line(11) = {3, 4};
Line(12) = {8, 7};

Curve Loop(1) = {7, 4, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {5, -11, 10, -1};
Plane Surface(2) = {2};
Curve Loop(3) = {10, 2, 8, 9};
Plane Surface(3) = {3};
Curve Loop(4) = {3, -7, 12, -8};
Plane Surface(4) = {4};
Curve Loop(5) = {4, 1, 2, 3};
Plane Surface(5) = {5};
Curve Loop(6) = {6, 12, 9, 11};
Plane Surface(6) = {6};

Surface Loop(1) = {5, 1, 4, 6, 3, 2};
Volume(1) = {1};

Physical Volume(1) = {1};
Physical Surface(2) = {1};
Physical Surface(3) = {3};
Physical Surface(4) = {1,2,3,4,5,6};
