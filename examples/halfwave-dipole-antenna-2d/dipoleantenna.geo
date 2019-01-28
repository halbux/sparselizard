length 		= 0.075;
width 		= 0.005;
radius		= 0.5;
gap 		= 0.005;

tsmall = 0.003;
tbig = 0.05;

Point(1) = {-width/2, gap/2, 0, tsmall};
Point(2) = {width/2, gap/2, 0, tsmall};
Point(3) = {width/2, gap/2+length, 0, tsmall};
Point(4) = {-width/2, gap/2+length, 0, tsmall};
Point(5) = {-width/2, -gap/2-length, 0, tsmall};
Point(6) = {width/2, -gap/2-length, 0, tsmall};
Point(7) = {width/2, -gap/2, 0, tsmall};
Point(8) = {-width/2, -gap/2, 0, tsmall};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {8, 1};
Line(10) = {7, 2};

Line Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
Line Loop(2) = {8, 5, 6, 7};
Plane Surface(2) = {2};
Line Loop(3) = {9, 1, -10, 7};
Plane Surface(3) = {3};

Point(9) = {0, 0, 0, tsmall};
Point(10) = {radius, 0, 0, tbig};
Point(11) = {-radius, 0, 0, tbig};
Point(12) = {0, radius, 0, tbig};
Point(13) = {0, -radius, 0, tbig};

Circle(11) = {11, 9, 13};
Circle(12) = {13, 9, 10};
Circle(13) = {10, 9, 12};
Circle(14) = {12, 9, 11};

Line Loop(4) = {14, 11, 12, 13};
Line Loop(5) = {8, 5, 6, 10, 2, 3, 4, -9};
Plane Surface(4) = {4, 5};


air = 1; conductor = 2; feed = 3; boundary = 4;

Physical Surface(air) = {4};
Physical Surface(conductor) = {1, 2};
Physical Surface(feed) = {3};
Physical Line(boundary) = {11, 12, 13, 14};
