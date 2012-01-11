// Gmsh project created on Thu Dec  2 14:47:53 2010
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {-1, 0, 0, 1.0};

Circle(1) = {4, 1, 3};
Circle(2) = {3, 1, 2};
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{2, 1};
}
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{3, 6};
}
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{12, 9};
}
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Line{15, 18};
}
Surface Loop(27) = {5, 11, 20, 26, 23, 17, 14, 8};
Volume(28) = {27};


Physical Surface(1) = {5, 26, 20, 11, 14, 8, 17, 23};
Physical Volume(2) = {28};
