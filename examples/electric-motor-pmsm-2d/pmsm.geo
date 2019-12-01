pi = 3.14159265359;

// Rotor [1] or stator [0] part of the geometry:
isrotor = 1;

// Central shaft radius [m]:
rshaft = 0.025;
// Rotor magnetic material radius [m]:
rrotmagmat = 0.055;
// Rotor permanent magnet radial length [m]:
drpermmag = 0.003;
// Permanent magnet angular coverage [degrees]:
alphapermmag = 33;

// Rotor-stator gap radial length [m]:
drrotstatgap = 0.001;

// Stator radius [m]:
rstat = 0.1;

// Stator winding slots (3 in total, one per phase):
wdopen = 0.004;
wdopendepth = 0.0012;
wddepth = 0.02;
wdwidthsmall = 0.006;
wdwidthlarge = 0.012;
// Curvature radius at the end of the winding slot:
rcurv = 0.003;

// Mesh sizes at shaft, permanent magnet, rotor-stator gap and stator:
msshaft = 0.006;
mspm = 0.0015;
msrsg = 0.0008;
msstat = 0.006;


// R is the current drawing radius:
R = 0;
Point(1) = {0, 0, 0, msshaft};

R = rshaft;
Point(2) = {R, 0, 0, msshaft};
Point(3) = {R*Cos(pi/4), R*Sin(pi/4), 0, msshaft};

R = rrotmagmat;
Point(6) = {R, 0, 0, mspm};
Point(7) = {R*Cos(pi/4), R*Sin(pi/4), 0, mspm};
Point(8) = {R*Cos( (pi/4+alphapermmag*pi/180)/2 ), R*Sin( (pi/4+alphapermmag*pi/180)/2 ), 0, mspm};
Point(9) = {R*Cos( (pi/4-alphapermmag*pi/180)/2 ), R*Sin( (pi/4-alphapermmag*pi/180)/2 ), 0, mspm};

R = rrotmagmat+drpermmag;
Point(10) = {R, 0, 0, msrsg};
Point(11) = {R*Cos(pi/4), R*Sin(pi/4), 0, msrsg};
Point(12) = {R*Cos( (pi/4+alphapermmag*pi/180)/2 ), R*Sin( (pi/4+alphapermmag*pi/180)/2 ), 0, msrsg};
Point(13) = {R*Cos( (pi/4-alphapermmag*pi/180)/2 ), R*Sin( (pi/4-alphapermmag*pi/180)/2 ), 0, msrsg};

R = rrotmagmat+drpermmag+drrotstatgap/2;
Point(14) = {R, 0, 0, msrsg};
Point(15) = {R*Cos(pi/4), R*Sin(pi/4), 0, msrsg};

R = rrotmagmat+drpermmag+drrotstatgap;
Point(16) = {R, 0, 0, msrsg};
Point(17) = {R*Cos(pi/4), R*Sin(pi/4), 0, msrsg};

R = rstat;
Point(18) = {R, 0, 0, msstat};
Point(19) = {R*Cos(pi/4), R*Sin(pi/4), 0, msstat};


// Stator winding slots (3 in total, one per phase):
For i In {0:2}

    alpharot = pi/4 / 3 /2 + i * pi/4 / 3;

    // Slot open:
    R = rrotmagmat+drpermmag+drrotstatgap;
    alpha1 = Atan(0.5*wdopen / R) + alpharot;
    alpha2 = -Atan(0.5*wdopen / R) + alpharot;
    Point(20+12*i) = {R*Cos(alpha1),R*Sin(alpha1),0, msrsg};
    Point(21+12*i) = {R*Cos(alpha2),R*Sin(alpha2),0, msrsg};

    R = rrotmagmat+drpermmag+drrotstatgap+wdopendepth;
    Point(22+12*i) = {R*Cos(alpha1),R*Sin(alpha1),0, msrsg};
    Point(23+12*i) = {R*Cos(alpha2),R*Sin(alpha2),0, msrsg};

    // Winding open:
    R = rrotmagmat+drpermmag+drrotstatgap+wdopendepth+0.5*wdopendepth;
    alpha1 = Atan(0.5*wdwidthsmall / R) + alpharot;
    alpha2 = -Atan(0.5*wdwidthsmall / R) + alpharot;
    Point(24+12*i) = {R*Cos(alpha1),R*Sin(alpha1),0, msrsg};
    Point(25+12*i) = {R*Cos(alpha2),R*Sin(alpha2),0, msrsg};

    R = rrotmagmat+drpermmag+drrotstatgap+wdopendepth+wddepth;
    alpha1 = Atan(0.5*wdwidthlarge / R) + alpharot;
    alpha2 = -Atan(0.5*wdwidthlarge / R) + alpharot;
    Point(26+12*i) = {R*Cos(alpha1),R*Sin(alpha1),0, msstat};
    Point(27+12*i) = {R*Cos(alpha2),R*Sin(alpha2),0, msstat};
    
    R = rrotmagmat+drpermmag+drrotstatgap+wdopendepth+wddepth;
    alpha1 = Atan((0.5*wdwidthlarge-rcurv) / R) + alpharot;
    alpha2 = -Atan((0.5*wdwidthlarge-rcurv) / R) + alpharot;
    Point(28+12*i) = {R*Cos(alpha1),R*Sin(alpha1),0, msstat};
    Point(29+12*i) = {R*Cos(alpha2),R*Sin(alpha2),0, msstat};
    
    R = rrotmagmat+drpermmag+drrotstatgap+wdopendepth+wddepth+rcurv;
    alpha1 = Atan((0.5*wdwidthlarge-rcurv) / R) + alpharot;
    alpha2 = -Atan((0.5*wdwidthlarge-rcurv) / R) + alpharot;
    Point(30+12*i) = {R*Cos(alpha1),R*Sin(alpha1),0, msstat};
    Point(31+12*i) = {R*Cos(alpha2),R*Sin(alpha2),0, msstat};

EndFor


Circle(1) = {2, 1, 3};
Circle(2) = {6, 1, 9};
Circle(3) = {9, 1, 8};
Circle(4) = {8, 1, 7};
Circle(5) = {10, 1, 13};
Circle(6) = {13, 1, 12};
Circle(7) = {12, 1, 11};
Circle(8) = {14, 1, 15};
Circle(9) = {16, 1, 21};
Circle(10) = {21, 1, 20};
Circle(11) = {20, 1, 33};
Circle(12) = {33, 1, 32};
Circle(13) = {32, 1, 45};
Circle(14) = {45, 1, 44};
Circle(15) = {44, 1, 17};
Circle(16) = {18, 1, 19};

Line(17) = {2, 6};
Line(18) = {6, 10};
Line(19) = {10, 14};
Line(20) = {14, 16};
Line(21) = {16, 18};
Line(22) = {3, 7};
Line(23) = {7, 11};
Line(24) = {11, 15};
Line(25) = {15, 17};
Line(26) = {17, 19};
Line(27) = {9, 13};
Line(28) = {8, 12};
Line(29) = {21, 23};
Line(30) = {23, 25};
Line(31) = {25, 27};
Line(32) = {31, 30};
Line(33) = {26, 24};
Line(34) = {24, 22};
Line(35) = {22, 20};
Line(36) = {33, 35};
Line(37) = {35, 37};
Line(38) = {37, 39};
Line(39) = {43, 42};
Line(40) = {38, 36};
Line(41) = {36, 34};
Line(42) = {34, 32};
Line(43) = {45, 47};
Line(44) = {47, 49};
Line(45) = {49, 51};
Line(46) = {55, 54};
Line(47) = {50, 48};
Line(48) = {48, 46};
Line(49) = {46, 44};
Line(50) = {23, 22};
Line(51) = {35, 34};
Line(52) = {47, 46};

Circle(53) = {27, 29, 31};
Circle(54) = {30, 28, 26};
Circle(55) = {39, 41, 43};
Circle(56) = {42, 40, 38};
Circle(57) = {51, 53, 55};
Circle(58) = {54, 52, 50};

Line Loop(1) = {1, 22, -4, -3, -2, -17};
Plane Surface(1) = {1};
Line Loop(2) = {2, 27, -5, -18};
Plane Surface(2) = {2};
Line Loop(3) = {3, 28, -6, -27};
Plane Surface(3) = {3};
Line Loop(4) = {4, 23, -7, -28};
Plane Surface(4) = {4};
Line Loop(5) = {5, 6, 7, 24, -8, -19};
Plane Surface(5) = {5};
Line Loop(6) = {8, 25, -15, -14, -13, -12, -11, -10, -9, -20};
Plane Surface(6) = {6};
Line Loop(7) = {10, -35, -50, -29};
Plane Surface(7) = {7};
Line Loop(8) = {12, -42, -51, -36};
Plane Surface(8) = {8};
Line Loop(9) = {14, -49, -52, -43};
Plane Surface(9) = {9};
Line Loop(10) = {50, -34, -33, -54, -32, -53, -31, -30};
Plane Surface(10) = {10};
Line Loop(11) = {51, -41, -40, -56, -39, -55, -38, -37};
Plane Surface(11) = {11};
Line Loop(12) = {52, -48, -47, -58, -46, -57, -45, -44};
Plane Surface(12) = {12};
Line Loop(13) = {21, 16, -26, -15, -49, -48, -47, -58, -46, -57, -45, -44, -43, -13, -42, -41, -40, -56, -39, -55, -38, -37, -36, -11, -35, -34, -33, -54, -32, -53, -31, -30, -29, -9};
Plane Surface(13) = {13};


// Physical region tags:
rotmagmat = 1; magnet = 2; magnetgap = 3; gaprot = 4; gapstat = 5; statmagmat = 6; windslot = 7; winda = 8; windb = 9; windc = 10;
gammarot = 11; gammastat = 12; gamma1rot = 13; gamma2rot = 14; gamma1stat = 15; gamma2stat = 16; inarc = 17; outarc = 18;

If (isrotor == 1)
    Physical Surface(1) = {1};
    Physical Surface(2) = {3};
    Physical Surface(3) = {2,4};
    Physical Surface(4) = {5};
    
    Physical Line(11) = {8};
    Physical Line(13) = {17,18,19};
    Physical Line(14) = {22,23,24};
    Physical Line(17) = {1};
Else
    Physical Surface(5) = {6};
    Physical Surface(6) = {13};
    Physical Surface(7) = {7,8,9};
    Physical Surface(8) = {10};
    Physical Surface(9) = {11};
    Physical Surface(10) = {12};
    
    Physical Line(12) = {8};
    Physical Line(15) = {20,21};
    Physical Line(16) = {25,26};
    Physical Line(18) = {16};
EndIf

