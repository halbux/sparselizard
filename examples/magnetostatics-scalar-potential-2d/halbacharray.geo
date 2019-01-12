// This file creates the 2D geometry of 5 aligned permanent magnets.

// Domain radius in meters:
rdom = 0.07;
// Magnet size in meters:
l = 0.01;
// Steel disk radius in meters:
rdisk = 0.01;

// Choose mesh size:
Mesh.CharacteristicLengthFactor = 0.1;


SetFactory("OpenCASCADE");
Rectangle(1) = {-2.5*l, 0, 0, l, l, 0};
Rectangle(2) = {-1.5*l, 0, 0, l, l, 0};
Rectangle(3) = {-0.5*l, 0, 0, l, l, 0};
Rectangle(4) = {0.5*l, 0, 0, l, l, 0};
Rectangle(5) = {1.5*l, 0, 0, l, l, 0};
Disk(6) = {3.5*l, rdisk+l*1.5, 0, rdisk, rdisk};

Disk(7) = {0, rdisk, 0, rdom, rdom};

// Break down the geometry:
Coherence;

// Define the physical regions:
magnet1 = 1; magnet2 = 2; magnet3 = 3; magnet4 = 4; magnet5 = 5; steel = 6; air = 7; zeropotential = 8;

Physical Surface(magnet1) = {1};
Physical Surface(magnet2) = {2};
Physical Surface(magnet3) = {3};
Physical Surface(magnet4) = {4};
Physical Surface(magnet5) = {5};
Physical Surface(steel) = {6};
Physical Surface(air) = {7};
Physical Point(zeropotential) = {2};


