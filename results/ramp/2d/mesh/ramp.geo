// Gmsh project created on Mon Sep 19 12:49:38 2011
R = 1.0;
lc = R/5;
Point(1) = {0, 0, 0, lc};
Point(2) = {-R, 0, 0, lc};
Point(3) = {R, 0, 0, lc};
Point(4) = {0, R, 0, lc};
Line(1) = {3, 1};
Line(2) = {1, 2};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 3};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};

Rotate {{0, 0, 1}, {0, 0, 0}, -0*Pi/4} {
  Surface{6};
}

Physical Line(2) = {4, 3};
Physical Line(3) = {1, 2};
Physical Surface(100) = {6};
Physical Point(20) = {2};
Physical Point(21) = {3};
