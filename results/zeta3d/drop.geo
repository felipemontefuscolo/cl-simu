// gmsh

lc = .2;
R = 1;

Point(1) = {0, 0, 0, lc};
Point(2) = {R, 0, 0, lc};
Point(3) = {0, R, 0, lc};
Point(4) = {0, 0, R, lc};
Point(5) = {0, 0,-R, lc};

Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {1, 4};

Line(4) = {1, 5};
Circle(5) = {2, 1, 3};
Circle(6) = {5, 1, 3};
Circle(7) = {4, 1, 3};
Circle(8) = {2, 1, 5};
Circle(9) = {2, 1, 4};


Line Loop(10) = -{7, -5, 9};
Ruled Surface(11) = {10};
Line Loop(12) = -{5, -6, -8};
Ruled Surface(13) = {12};


Line Loop(14) = -{2, -7, -3};
Plane Surface(15) = {14};
Line Loop(16) = {2, -6, -4};
Plane Surface(17) = {16};
Line Loop(18) = -{1, 3, -9};
Plane Surface(19) = {18};
Line Loop(20) = {1, 4, -8};
Plane Surface(21) = {20};



//Transfinite Line {3, 4, 1, 2, 5, 6, 7, 8, 9} = 2 Using Progression 1;

Physical Line(1) = {3, 4}; // corner
Physical Line(2) = {7, 6}; // wall
Physical Line(5) = {9, 8}; // floor

Physical Surface(2) = {15, 17}; // wall
Physical Surface(4) = {19, 21}; // floor
Physical Surface(3) = {13, 11}; // interface


Physical Point(1) = {4, 5};


//Line Loop(22) = {5, -2, -1};
//Plane Surface(23) = {22};
//Surface Loop(24) = {23, 15, 19, 11};
//Volume(25) = {24};
//Surface Loop(26) = {13, 21, 17, 23};
//Volume(27) = {26};
//Physical Volume(123) = {25, 27};

Surface Loop(22) = {13, 11, 15, 17, 21, 19};
Volume(23) = {22};
Physical Volume(123) = {23};


// use this to force a coarse mesh inside
//Mesh.CharacteristicLengthExtendFromBoundary = 1;
//Mesh.CharacteristicLengthMax = 0.2;



