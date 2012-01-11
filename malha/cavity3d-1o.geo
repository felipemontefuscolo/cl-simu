cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Point(5) = {1, 0, 1, cl1};
Point(6) = {1, 1, 1, cl1};
Point(10) = {0, 1, 1, cl1};
Point(14) = {0, 0, 1, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(8) = {5, 6};
Line(9) = {6, 10};
Line(10) = {10, 14};
Line(11) = {14, 5};
Line(13) = {2, 5};
Line(14) = {3, 6};
Line(18) = {4, 10};
Line(22) = {1, 14};
Line Loop(6) = {2, 3, 4, 1};
Plane Surface(6) = {6};
Line Loop(15) = {2, 14, -8, -13};
Ruled Surface(15) = {15};
Line Loop(19) = {3, 18, -9, -14};
Ruled Surface(19) = {19};
Line Loop(23) = {4, 22, -10, -18};
Ruled Surface(23) = {23};
Line Loop(27) = {1, 13, -11, -22};
Ruled Surface(27) = {27};
Line Loop(28) = {8, 9, 10, 11};
Plane Surface(28) = {28};
Surface Loop(1) = {6, 28, 15, 19, 23, 27};
Volume(1) = {1};

// condições de contorno
Physical Point(2) = {3, 4, 6, 10};
Physical Line(2) = {3, 9, 14, 18, 2, 4, 8, 10, 1, 22, 11, 13};
Physical Surface(1) = {19};
Physical Surface(2) = {6, 15, 23, 27, 28};
Physical Volume(3) = {1};


Mesh.Smoothing = 100;


// TRANSFINITE
Transfinite Line {4, 10, 2, 8, 22, 13, 14, 18, 3, 9, 11, 1} = 17 Using Progression 1;
Transfinite Surface {23} = {4, 1, 14, 10};
Transfinite Surface {15} = {3, 2, 5, 6};
Transfinite Surface {28} = {10, 14, 5, 6};
Transfinite Surface {6} = {4, 1, 2, 3};
Transfinite Surface {19} = {4, 10, 6, 3};
Transfinite Surface {27} = {1, 14, 5, 2};

Transfinite Volume{1} = {4, 1, 14, 10, 3, 2, 5, 6};
