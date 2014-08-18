lc = .4;
R = 1;




Point(1) = {0, 0, 0, lc};
Point(2) = {R, 0, 0, lc};
Point(3) = {0, R, 0, lc};
Point(4) = {0, 0, R, lc};
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {1, 4};
Circle(4) = {2, 1, 3};

Circle(5) = {4, 1, 3};
Circle(6) = {4, 1, 2};
Line Loop(7) = {4, -2, -1};
Plane Surface(8) = {7};
Line Loop(9) = {5, -2, 3};
Plane Surface(10) = {9};
Line Loop(11) = {1, 3, 6};
Plane Surface(12) = {11};


Line Loop(13) = {5, -4, -6};
Ruled Surface(14) = {13} In Sphere {1};



Surface Loop(15) = {14, 10, 8, 12};
Volume(16) = {15};

Physical Point(1) = {1}; // origin
Physical Line(2) = {6};  // cl
Physical Surface(3) = {14}; // interf
Physical Surface(4) = {10}; // solid
Physical Surface(5) = {12}; // solid
Physical Surface(6) = {8};  // solid
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Point(7) = {2};
Physical Point(8) = {3};
Physical Point(9) = {4};

Transfinite Line {2, 3, 1, 6, 5, 4} = 2 Using Progression 1;

Physical Volume(123) = {16};


Transfinite Surface {14} = {3, 4, 2} Left;
Transfinite Surface {10} = {3, 1, 4} Left;
Transfinite Surface {8} = {3, 2, 1} Left;
Transfinite Surface {12} = {1, 2, 4} Left;



Physical Line(4) = {5};
Physical Line(6) = {4};




