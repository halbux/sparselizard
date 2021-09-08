Mesh.CharacteristicLengthFactor = 0.15;

SetFactory("OpenCASCADE");

lx = 7;
ly = 10;

Rectangle(1) = {1e-8, -ly/2, 0, lx, ly, 0};
Rectangle(2) = {1, -2.25, 0, 0.5, 4.5, 0};
Disk(3) = {3, 0.5, 0, 0.5, 0.5};
Disk(4) = {1.25, 2.25, 0, 0.25, 0.25};
Disk(5) = {1.25, -2.25, 0, 0.25, 0.25};

Coherence;

Physical Surface(1) = {5,6,7,8,9};
Physical Surface(2) = {3};
Physical Surface(3) = {4};
Physical Line(4) = {8};
Physical Line(5) = {1,3,4};
