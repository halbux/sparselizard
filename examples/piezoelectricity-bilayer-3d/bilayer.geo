// PZT thickness [um]:
pztthickness = 10;
// Silicon thickness [um]:
siliconthickness = 20;

// Bilayer length [um]:
length = 200;
// Bilayer width [um]:
width = 100;


SetFactory("OpenCASCADE");
Box(1) = {0, 0, siliconthickness, length, width, pztthickness};
Box(2) = {0, 0, 0, length, width, siliconthickness};

Physical Volume(1) = 1;
Physical Volume(2) = 2;

Physical Surface(3) = 6;
Physical Surface(4) = 12;
Physical Surface(5) = 11;

Physical Surface(6) = {7,1};
Physical Surface(7) = {8,2};

Physical Surface(8) = {1,2,3,4,5,6};
Physical Surface(9) = {7,8,9,10,11,1,2,3,4,6};

// Mesh scaling factor:
Mesh.ScalingFactor = 1e-6;
