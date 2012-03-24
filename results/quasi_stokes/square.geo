// square domain

cl1 = 1;
cl2 = 1;

Point(1) = {  0,   0, 0, cl2};
Point(2) = {  0, 0.5, 0, cl2};
Point(3) = {  0,   1, 0, cl1};
Point(4) = {0.5,   1, 0, cl1};
Point(5) = {  1,   1, 0, cl1};
Point(6) = {  1, 0.5, 0, cl2};
Point(7) = {  1,   0, 0, cl2};
Point(8) = {0.5,   0, 0, cl2};

Point(9) = {0.5, 0.5, 0, cl2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Line(9) = {4, 9};
Line(10) = {9, 2};
Line(11) = {9, 8};
Line(12) = {9, 6};

Line Loop(13) = {2, 3, 9, 10};
Plane Surface(14) = {13};
Line Loop(15) = {4, 5, -12, -9};
Plane Surface(16) = {15};
Line Loop(17) = {11, 8, 1, -10};
Plane Surface(18) = {17};
Line Loop(19) = {12, 6, 7, -11};
Plane Surface(20) = {19};

Physical Point(2) = {3, 5};
Physical Line(2) = {2, 1, 8, 7, 6, 5};
Physical Line(1) = {3, 4};
Physical Line(3) = {9, 10, 11, 12};
Physical Surface(3) = {14, 18, 20, 16};

//Mesh.Smoothing = 100;


Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12} = 2 Using Progression 1;

Transfinite Surface {18} Right;
Transfinite Surface {16} Left;
Transfinite Surface {14} Left;
Transfinite Surface {20} Right;


