// Gmsh project created on Mon Sep 19 12:49:38 2011
R = 1.0;
lc = R/5;
Point(1) = {0, 0, 0, lc};
Point(2) = {-R, 0, 0, 0.25*lc};
Point(3) = {R, 0, 0, 0.25*lc};
Point(4) = {0, R, 0, lc};


Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/2.25} {
  Duplicata { Point{3}; }
}
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2.25} {
  Duplicata { Point{2}; }
}
Circle(1) = {6, 1, 2};
Circle(2) = {2, 1, 4};
Circle(3) = {4, 1, 3};
Circle(4) = {3, 1, 5};
Line(5) = {5, 6};
Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(7) = {6};

Physical Point(20) = {6};
Physical Point(21) = {5};
Physical Line(2) = {1, 2, 3, 4};
Physical Line(3) = {5};
Physical Surface(100) = {7};
