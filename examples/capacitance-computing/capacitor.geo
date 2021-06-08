// Characteristic mesh size:
tcfine = 0.2e-3;
tccoarse = 2e-3;

hdielec = 2e-3;
welectrode = 5e-3;
wbox = 40e-3;
hbox = 20e-3;

Point(1) = {-welectrode/2, 0, 0, tcfine};
Point(2) = {welectrode/2, 0, 0, tcfine};
Point(3) = {-wbox/2, 0, 0, tccoarse};
Point(4) = {wbox/2, 0, 0, tccoarse};

Point(5) = {-welectrode/2, hdielec, 0, tcfine};
Point(6) = {welectrode/2, hdielec, 0, tcfine};
Point(7) = {-wbox/2, hdielec, 0, tccoarse};
Point(8) = {wbox/2, hdielec, 0, tccoarse};

Point(9) = {-wbox/2, hbox, 0, tccoarse};
Point(10) = {wbox/2, hbox, 0, tccoarse};

Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {2, 4};
Line(4) = {4, 8};
Line(5) = {8, 10};
Line(6) = {10, 9};
Line(7) = {9, 7};
Line(8) = {7, 3};
Line(9) = {7, 5};
Line(10) = {5, 6};
Line(11) = {6, 8};

Line Loop(1) = {8, 1, 2, 3, 4, -11, -10, -9};
Plane Surface(1) = {1};
Line Loop(2) = {7, 9, 10, 11, 5, 6};
Plane Surface(2) = {2};

dielectric = 1; air = 2; electrode = 3; ground = 4;

Physical Surface(dielectric) = 1;
Physical Surface(air) = 2;
Physical Line(electrode) = 10;
Physical Line(ground) = {1,2,3};
