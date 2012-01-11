// Gmsh project created on Thu Dec  2 20:26:13 2010
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6};
}
Physical Surface(1) = {19};
Physical Surface(2) = {6, 23, 28, 27, 15};
Physical Line(2) = {9, 14, 3, 18};
Physical Point(2) = {6, 3, 4, 10};
Physical Volume(3) = {1};
