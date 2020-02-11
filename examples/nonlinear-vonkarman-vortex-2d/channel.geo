boxlen = 2.0;
boxwidth = 0.8;

cylinderradius = 0.07;
cylinderxpos = 0.25;
// Tiny misalignment:
cylinderypos = boxwidth/2 + boxwidth/1000;

lc = 0.02;
lcfine = 0.002;

Point(1) = {0,0,0, lc};
Point(2) = {boxlen,0,0, lc};
Point(3) = {boxlen,boxwidth,0, lc};
Point(4) = {0,boxwidth,0, lc};

Point(5) = {cylinderxpos,cylinderypos,0, lcfine};
Point(6) = {cylinderxpos+cylinderradius,cylinderypos,0, lcfine};
Point(7) = {cylinderxpos,cylinderypos+cylinderradius,0, lcfine};
Point(8) = {cylinderxpos-cylinderradius,cylinderypos,0, lcfine};
Point(9) = {cylinderxpos,cylinderypos-cylinderradius,0, lcfine};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Line Loop(1) = {4, 1, 2, 3};
Line Loop(2) = {6, 7, 8, 5};
Plane Surface(1) = {1, 2};
Plane Surface(2) = {2};


fluid = 1; cylinder = 2; inlet = 3; outlet = 4;

Physical Surface(fluid) = 1;
Physical Surface(cylinder) = 2;
Physical Line(inlet) = 4;
Physical Line(outlet) = 2;
