// Gmsh project created on Fri May 27 23:13:19 2011
// [-L, 1.0] x [-L, L]
//
//
//

d = 10;

L = 1.;

Point(1) = { 0,  0, 0, d};
Point(2) = { L,  0, 0, d};
Point(3) = { L,  L+.1, 0, d};
Point(4) = { 0,  L, 0, d};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


//Reverse Surface { 10,16,14,12 };


Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Line(2) = {2, 3, 4};
Physical Line(1) = {1};
Physical Surface(92) = {6};
