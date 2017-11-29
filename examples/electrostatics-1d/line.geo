// Choose the mesh size:
c = 0.1;

Point(1) = {0, 0, 0, c};
Point(2) = {1, 0, 0, c};

Line(1) = {1, 2};

Physical Line(1) = {1};
Physical Point(2) = {1};
Physical Point(3) = {2};
