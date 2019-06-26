Mesh.ScalingFactor = 1e-6;

tcfine = 1;

Mesh.CharacteristicLengthMin = 0.5;
Mesh.CharacteristicLengthMax = 4;

SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 320, 120, 0};

Point(5) = {100, 0, 0, 100};
Point(6) = {100, 60, 0, tcfine};
Point(7) = {102, 60, 0, tcfine};
Point(8) = {102, 62, 0, tcfine};
Point(9) = {104, 60, 0, tcfine};
Point(10) = {104, 0, 0, 100};
Line(5) = {5, 6};
Line(6) = {9, 10};
Line(9) = {5, 10};

Circle(7) = {6, 7, 8};
Circle(8) = {8, 7, 9};

Line Loop(2) = {1, 2, 3, 4};
Line Loop(3) = {5, 7, 8, 6, -9};
Plane Surface(2) = {3};

Coherence;


fluid = 1; pillar = 2; inlet = 3; outlet = 4; sides = 5; clamp = 6;

Physical Surface(fluid) = {3};
Physical Surface(pillar) = {2};
Physical Line(inlet) = {1};
Physical Line(outlet) = {3};
Physical Line(sides) = {2,4,9};
Physical Line(clamp) = {10};
