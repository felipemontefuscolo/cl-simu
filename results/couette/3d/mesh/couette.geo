// Gmsh project created on Mon Feb 27 16:51:38 2012

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {0.0, 0.0, 0.5, 1.0};
Point(3) = {0.0, 0.5, 0.0, 1.0};
Point(4) = {0.0, 0.5, 0.5, 1.0};

Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Extrude {0.5, 0, 0} {
  Surface{6};
}


Physical Point(2) = {6, 10, 14, 5};
Physical Line(2) = {8, 9, 10, 11};
Physical Surface(2) = {27, 15, 6, 19, 23};
Physical Surface(1) = {28};
Physical Volume(33) = {1};