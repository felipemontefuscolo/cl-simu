// Gmsh project created on Wed Oct 19 14:36:43 2011
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {-1, 0, 0, 1.0};
Point(4) = {0, 0, -1, 1.0};


Point(5) = {0, 0, 1, 1.0};
Point(6) = {0, 1, 0, 1.0};



Circle(1) = {5, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 2};
Circle(4) = {2, 1, 5};
Circle(5) = {3, 1, 6};
Circle(6) = {4, 1, 6};
Circle(7) = {2, 1, 6};
Circle(8) = {5, 1, 6};
Line Loop(9) = {2, 6, -5};
Ruled Surface(10) = {9};
Line Loop(11) = {6, -7, -3};
Ruled Surface(12) = {11};
Line Loop(13) = {7, -8, -4};
Ruled Surface(14) = {13};
Line Loop(15) = {8, -5, -1};
Ruled Surface(16) = {15};
Line Loop(17) = {2, 3, 4, 1};
Plane Surface(18) = {17};
Surface Loop(19) = {18, 10, 12, 14, 16};
Volume(20) = {19};

// linha triplice
Physical Line(20) = {1, 4, 3, 2};
// solido-fluido
Physical Surface(3) = {18};
// fluid0-gas
Physical Surface(2) = {14, 12, 10, 16};
// fluido interior
Physical Volume(100) = {20};








