// Gmsh project created on Mon Oct 25 15:56:52 2010

N = 20;
R = 1;
A = 0.2;
n = 2;   // modo de oscilacao
lc = R/5;


d = 1;

// tau = R^2 
Point(1) = {0, 0, 0, d};
Point(2) = {R, 0, 0, d};
Point(3) = {-R, 0, 0, d};
Point(4) = {0, R, 0, d};
Point(5) = {0, -R, 0, d};


Circle(1) = {5, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 2};
Circle(4) = {2, 1, 5};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Physical Line(2) = {3, 2, 1, 4};
Physical Surface(88) = {6};
Physical Point(2) = {2, 5, 3, 4};
