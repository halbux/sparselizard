tc = 0.1;

Point(1) = {0, 0, 0, tc};
Point(2) = {1, 0, 0, tc};
Point(3) = {0, 1, 0, tc};
Point(4) = {1, 1, 0, tc};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

// Try to get a tri/quad mix in the mesh:
Transfinite Line {1} = 1/tc;

Line Loop(5) = {4, 1, 2, 3};
Plane Surface(1) = {5};

Recombine Surface {1};


Physical Surface(1) = {1};
Physical Line(2) = {4};
Physical Line(3) = {2};
Physical Line(4) = {1,2,3,4};
