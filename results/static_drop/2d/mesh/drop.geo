// Gmsh project created on Mon Oct 25 15:56:52 2010

// use smoothing = 100

R = 1;
xc = 0;
yc = 0;

lc = R/3;


Point(1) = {xc, yc, 0, lc};
Point(2) = {xc-R, yc, 0, lc};
Point(3) = {xc, yc+R, 0, lc};
Point(4) = {xc+R, yc, 0, lc};
Point(5) = {xc, yc-R, 0, lc};

Circle(1) = {5, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 4};
Circle(4) = {4, 1, 5};
Line Loop(5) = {3, 4, 1, 2};

Line(7) = {3, 1};
Line(8) = {1, 5};


Line Loop(9) = {7, 8, -4, -3};
Plane Surface(10) = {-9};
Line Loop(11) = {7, 8, 1, 2};
Plane Surface(12) = {11};

i = 2; // 0 1 2 3 4 5 ...
n0 = 5;
n = (2^i)*(n0-1)+1;
Transfinite Line {3, 4, 1, 2} = n Using Bump 0.2;
Transfinite Line {8, 7} = n Using Bump 1;
Transfinite Surface {10} Left;
Transfinite Surface {12} Left;




Physical Line(2) = {2, 3, 1, 4};
Physical Line(7) = {7, 8};
Physical Surface(5) = {12, 10};



