// Waveguide width and height [um]:
w = 0.5;
h = 0.25;
// Distance between waveguides [um]:
d = 0.1;

// Outer box width and height:
bw = 4;
bh = 2;

tw = 0.02;
tout = 0.2;


Point(1) = {-d/2-w, -h/2, 0, tw};
Point(2) = {-d/2, -h/2, 0, tw};
Point(3) = {-d/2, h/2, 0, tw};
Point(4) = {-d/2-w, h/2, 0, tw};

Point(5) = {d/2, -h/2, 0, tw};
Point(6) = {d/2+w, -h/2, 0, tw};
Point(7) = {d/2+w, h/2, 0, tw};
Point(8) = {d/2, h/2, 0, tw};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Point(9) = {-bw/2, -bh/2, 0, tout};
Point(10) = {bw/2, -bh/2, 0, tout};
Point(11) = {bw/2, bh/2, 0, tout};
Point(12) = {-bw/2, bh/2, 0, tout};

Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

Line Loop(3) = {12, 9, 10, 11};
Line Loop(4) = {4, 1, 2, 3};
Line Loop(5) = {5, 6, 7, 8};
Plane Surface(3) = {3, 4, 5};


// Mesh scaling factor:
Mesh.ScalingFactor = 1e-6;

// Left waveguide:
Physical Surface(1) = {1};
// Right waveguide:
Physical Surface(2) = {2};
// Clad:
Physical Surface(3) = {3};
// Whole domain:
Physical Surface(4) = {1,2,3};



