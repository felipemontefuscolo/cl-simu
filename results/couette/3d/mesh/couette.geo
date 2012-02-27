// Gmsh project created on Mon Feb 27 16:51:38 2012

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {0.0, 0.0, 0.5, 1.0};
Point(3) = {0.0, 0.5, 0.0, 1.0};
Point(4) = {0.0, 0.5, 0.5, 1.0};
Point(5) = {0.5, 0.0, 0.0, 1.0};
Point(6) = {0.5, 0.0, 0.5, 1.0};
Point(7) = {0.5, 0.5, 0.0, 1.0};
Point(8) = {0.5, 0.5, 0.5, 1.0};


Line(1) = {2, 6};
Line(2) = {6, 5};
Line(3) = {5, 1};
Line(4) = {1, 2};
Line(5) = {4, 3};
Line(6) = {3, 7};
Line(7) = {7, 8};
Line(8) = {8, 4};
Line(9) = {7, 5};
Line(10) = {6, 8};
Line(11) = {1, 3};
Line(12) = {4, 2};
Line Loop(13) = {11, -5, 12, -4};
Plane Surface(14) = {13};
Line Loop(15) = {6, 7, 8, 5};
Plane Surface(16) = {15};
Line Loop(17) = {9, -2, 10, -7};
Plane Surface(18) = {17};
Line Loop(19) = {3, 4, 1, 2};
Plane Surface(20) = {19};
Line Loop(21) = {9, 3, 11, 6};
Plane Surface(22) = {21};
Line Loop(23) = {10, 8, 12, 1};
Plane Surface(24) = {23};
Surface Loop(25) = {24, 18, 22, 20, 14, 16};
Volume(26) = {25};




Physical Volume(27) = {26};
Physical Surface(1) = {18};
Physical Surface(2) = {22, 16, 14, 20, 24};
Physical Point(2) = {8, 6, 5, 7};
