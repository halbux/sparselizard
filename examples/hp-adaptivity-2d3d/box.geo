Mesh.CharacteristicLengthFactor = 5;

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};

Physical Volume(1) = {1};
Physical Surface(2) = {1};
Physical Surface(3) = {2};
