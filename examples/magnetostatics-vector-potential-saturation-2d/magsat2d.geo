// Choose the mesh size:
Mesh.CharacteristicLengthFactor = 8;
Mesh.MeshSizeFromCurvature = 400;

SetFactory("OpenCASCADE");
Disk(1) = {-0.0125, 0, 0, 0.01, 0.01};
Disk(2) = {-0.0125, 0, 0, 0.011, 0.011};
Disk(3) = {0.0125, 0, 0, 0.01, 0.01};
Disk(4) = {0.0125, 0, 0, 0.011, 0.011};
Disk(5) = {0, 0, 0, 0.03, 0.02};
Disk(6) = {0, 0, 0, 0.1, 0.1};

Coherence;

conductor1 = 1; conductor2 = 2; insulator = 3; steel = 4; air = 5;

Physical Surface(conductor1) = {1};
Physical Surface(conductor2) = {3};
Physical Surface(insulator) = {4,5};
Physical Surface(steel) = {6};
Physical Surface(air) = {7};
