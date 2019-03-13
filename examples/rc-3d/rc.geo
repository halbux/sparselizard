// This defines the geometry for a rectangular conductor trace connected to a circular 
// shaped capacitor made of a dielectric layer sandwiched between two conducting traces.

// Dimensions in um:
thicknessdielectric = 30;
thicknesstraces = 1;
diskradius = 300;
lengthtrace = diskradius+200;
widthtrace = 50;

// Set the characteristic mesh size required (a finer mesh gives more accurate results but takes longer to simulate):
Mesh.CharacteristicLengthFactor = 0.2;


SetFactory("OpenCASCADE");

// The bottom trace of the disk:
Cylinder(1) = {0, 0, 0, 0, 0, thicknesstraces, diskradius, 2*Pi};
// The dielectric of the disk:
Cylinder(2) = {0, 0, thicknesstraces, 0, 0, thicknessdielectric, diskradius, 2*Pi};
// The top trace of the disk:
Cylinder(3) = {0, 0, thicknesstraces+thicknessdielectric, 0, 0, thicknesstraces, diskradius, 2*Pi};
// The bottom trace leading to the disk:
Box(4) = {-lengthtrace, -widthtrace/2, 0, lengthtrace-diskradius/2, widthtrace, thicknesstraces};
// The top trace leading to the disk:
Box(5) = {-lengthtrace, -widthtrace/2, thicknesstraces+thicknessdielectric, lengthtrace-diskradius/2, widthtrace, thicknesstraces};

Coherence;

// All dimensions above are in um:
Mesh.ScalingFactor = 1e-6;


// Define the regions of interest for the simulation:
electrode = 1; ground = 2; conductor = 3; dielectric = 4; 

Physical Surface(electrode) = {26};
Physical Surface(ground) = {21};
Physical Volume(conductor) = {3,4,5,6,7,8};
Physical Volume(dielectric) = {2};
