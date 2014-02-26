Point(1) = {4, 1, 0, 1.0};
Point(2) = {4, 0, 0, 1.0};
Point(3) = {12, 0, 0, 1.0};
Point(4) = {12, 2, 0, 1.0};
Point(5) = {0, 2, 0, 1.0};
Point(6) = {0, 1, 0, 1.0};

Point(7) = {4, 2, 0, 1.0};


Line(1) = {6, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 7};
Line(6) = {7, 1};
Line(7) = {7, 5};
Line(8) = {5, 6};
Line Loop(9) = {5, 6, 2, 3, 4};
Plane Surface(10) = {9};
Line Loop(11) = {7, 8, 1, -6};
Plane Surface(12) = {11};

Transfinite Line {4} = 31 Using Progression 1;
Transfinite Line {2, 6, 8} = 16 Using Progression 1;
Transfinite Line {7, 1} = 61 Using Progression 1;
Transfinite Line {5, 3} = 121 Using Progression 1;

Transfinite Surface {10} = {7, 2, 3, 4};
Transfinite Surface {12} = {5, 6, 1, 7};


Physical Line(1) = {8};
Physical Line(2) = {4};
Physical Line(3) = {5, 7, 1, 3, 2};



Physical Surface(13) = {12, 10};
