// Increase n to have a finer structured mesh:
n = 10;


Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};

Transfinite Line {1} = n;
Transfinite Line {2} = n;
Transfinite Line {3} = n;
Transfinite Line {4} = n;

Line Loop(5) = {4, 1, 2, 3};
Plane Surface(1) = {5};

Transfinite Surface {1} = {1,2,3,4};
Recombine Surface {1};


Physical Surface(1) = {1};
