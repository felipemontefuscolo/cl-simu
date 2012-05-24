// Gmsh project created on Mon Sep 19 12:49:38 2011
R = 0.125;
lc = R/5;
Point(1) = {0, 0, 0, lc};
Point(2) = {-R, 0, 0, lc};
Point(3) = {R, 0, 0, lc};
Point(4) = {0, R, 0, lc};
Line(1) = {3, 1};
Line(2) = {2, 1};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 3};


hc = 2000;

Point(5) = {0, 0.7*R, 0, hc};
Point(6) = {0, 0.2*R, 0, hc};
Point(7) = {-0.6*R, 0.3*R, 0, hc};
Point(8) = {0.6*R, 0.3*R, 0, hc};


Spline(5) = {7, 5, 8, 6, 7};




Line Loop(6) = {3, 4, 1, -2};
Line Loop(7) = {5};
Plane Surface(8) = {6, 7};
Plane Surface(9) = {7};



Transfinite Line {2, 1} = 70/1 Using Progression 1.02; //baixo FINE MESH 750
Transfinite Line {4, 3} = 197/1 Using Progression 1; // cima
Transfinite Line {5} = 42/1 Using Progression 1; // meio




Physical Point(20) = {2};
Physical Point(21) = {3};
Physical Line(3) = {2, 1};
Physical Line(2) = {4, 3};
Physical Surface(100) = {8, 9};


//Physical Point(23) = {1};

