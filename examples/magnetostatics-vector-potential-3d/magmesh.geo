// Wire position and radius:
r = 0.1;
xpos = -0.5;
// Shield position and radius:
xposshield = 0.7;
rshieldin = 0.35;
rshieldout = 0.5;
// Domain size around:
rout = 1.5;
// Height:
h = 0.3;

SetFactory("OpenCASCADE");
Disk(1) = {xpos, 0, 0, r, r};
Disk(2) = {xposshield, 0, 0, rshieldout, rshieldout};
Disk(3) = {xposshield, 0, 0, rshieldin, rshieldin};
Disk(4) = {0, 0, 0, rout, rout};

Characteristic Length {4, 2, 1} = 0.04;
Characteristic Length {3} = 0.04;

Coherence;

Recombine Surface(1);
Recombine Surface(3);
Recombine Surface(4);
Recombine Surface(5);

nz = 1;
Extrude {0,0,h} { Surface{1}; Layers{nz}; Recombine;}
Extrude {0,0,h} { Surface{3}; Layers{nz}; Recombine;}
Extrude {0,0,h} { Surface{4}; Layers{nz}; Recombine;}
Extrude {0,0,h} { Surface{5}; Layers{nz}; Recombine;}


// Define the physical regions:
conductor = 1; shield = 2; air = 3; contour = 4;

Physical Volume(conductor) = {1};
Physical Volume(air) = {2,4};
Physical Volume(shield) = {3};
Physical Surface(contour) = {13, 1,3,4,5, 7,9,12,16};
