// cavity
cl1 = 0.1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Point(5) = {1, 0, 1, cl1};
Point(6) = {1, 1, 1, cl1};
Point(7) = {0, 1, 1, cl1};
Point(8) = {0, 0, 1, cl1};


Line(1) = {8, 7};
Line(2) = {7, 6};
Line(3) = {6, 5};
Line(4) = {5, 8};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {7, 4};
Line(10) = {3, 6};
Line(11) = {1, 8};
Line(12) = {5, 2};
Line Loop(13) = {4, 1, 2, 3};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};
Line Loop(17) = {11, 1, 9, 8};
Plane Surface(18) = {17};
Line Loop(19) = {12, 6, 10, 3};
Plane Surface(20) = {-19};
Line Loop(21) = {2, -10, 7, -9};
Plane Surface(22) = {-21};
Line Loop(23) = {4, -11, 5, -12};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 24, 18, 22, 20, 16};
Volume(26) = {25};

Physical Surface(1) = {22};
Physical Surface(2) = {14, 16, 18, 20, 24};
Physical Line(2) = {7, 9, 2, 10};
Physical Volume(9) = {26};

Mesh.Smoothing = 100;

Transfinite Line {7, 2, 9, 8, 1, 11, 5, 4, 3, 12, 6, 10} = 7 Using Progression 1;

Transfinite Surface {14} = {8, 7, 6, 5};
Transfinite Surface {16} = {1, 4, 3, 2};
Transfinite Surface {18} = {8, 7, 4, 1};
Transfinite Surface {20} = {5, 6, 3, 2};
Transfinite Surface {22} = {7, 6, 3, 4};
Transfinite Surface {24} = {8, 5, 2, 1};

Transfinite Volume{26} = {8, 7, 6, 5, 1, 4, 3, 2};






