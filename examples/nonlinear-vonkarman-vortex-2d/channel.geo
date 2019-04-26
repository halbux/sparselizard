boxlen = 2.0;
boxwidth = 0.8;

cylinderradius = 0.07;
cylinderxpos = 0.25;
// Tiny misalignment:
cylinderypos = boxwidth/2 + boxwidth/1000;

Mesh.CharacteristicLengthFactor = 0.05;


SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, boxlen, boxwidth, 0};
Disk(2) = {cylinderxpos, cylinderypos, 0, cylinderradius, cylinderradius};

Coherence;

fluid = 1; cylinder = 2; inlet = 3; outlet = 4;

Physical Surface(fluid) = 3;
Physical Surface(cylinder) = 2;
Physical Line(inlet) = 2;
Physical Line(outlet) = 3;
