// Gmsh project created on Fri May 27 23:13:19 2011
// [-L, 1.0] x [-L, L]
//
//
//

d = 1;

L = 1.2;

Point(1) = { 0,  0, 0, d};
Point(2) = {-L, -L, 0, d};
Point(3) = {-L, +L, 0, d};
Point(4) = {+L, +L, 0, d};
Point(5) = {+L, -L, 0, d};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};


// h = L/(N-1);

// meshes in paper:
// 19x13
// 31x21
// 43x29
// 61x41

Line(5) = {1, 3};
Line(6) = {1, 4};
Line(7) = {1, 5};
Line(8) = {1, 2};
Line Loop(9) = {5, -1, -8};
Plane Surface(10) = -{9};
Line Loop(11) = {6, -2, -5};
Plane Surface(12) = -{11};
Line Loop(13) = {3, -7, 6};
Plane Surface(14) = {13};
Line Loop(15) = {7, 4, -8};
Plane Surface(16) = {15};


Physical Surface(4) = {10, 16, 14, 12};
Physical Line(2) = {2, 1, 4};
Physical Line(1) = {3};
Physical Point(2) = {4, 5};

Reverse Surface { 10,16,14,12 };



