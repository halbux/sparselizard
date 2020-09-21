Mesh.CharacteristicLengthFactor = 1.0;

SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 1, 1, 0};
Disk(2) = {0.5, 0.15, 0, 0.02, 0.02};

Coherence;

fluid = 1;
disk = 2;
inlet = 3;
outlet = 4;
sides = 5;
    
Physical Surface(disk) = {2};
Physical Surface(fluid) = {3};

Physical Line(inlet) = {1};
Physical Line(outlet) = {4};
Physical Line(sides) = {2,3};
