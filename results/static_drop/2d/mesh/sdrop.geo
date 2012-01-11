// Gmsh project created on Mon Oct 25 15:56:52 2010

// use MeshAdapt algorithm

R = 1;
xc = 0;
yc = 0;

lc = R/3;

Point(1) = {xc, yc, 0, lc};
Point(2) = {xc-R, yc, 0, lc};
Point(3) = {xc, yc+R, 0, lc};
Point(4) = {xc+R, yc, 0, lc};
Point(5) = {xc, yc-R, 0, lc};

Circle(1) = {5, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 4};
Circle(4) = {4, 1, 5};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

Physical Line(2) = {3, 1, 4, 2};
Physical Surface(5) = {6};

//Mesh.Smoothing = 12;
