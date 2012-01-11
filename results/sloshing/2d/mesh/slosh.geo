// Gmsh project created on Fri Sep  2 15:21:47 2011

h = 1.5;
d = 1.0;
a = 0.1;

// mesh
lc = d/7;

Point(1) = {0, 0, 0, lc};
Point(2) = {d, 0, 0, lc};

Point(3) = {0, h+a, 0, lc};
Point(4) = {d/40, h+a, 0, lc};
N=10;
For s In {1:N}
	
	x = s*d/(N+1);
	y = h + a*Cos(Pi*x/d);
	
	Point(4+s) = {x,y,0,lc};

EndFor
Point(5+N) = {d*(1-1/40), h-a, 0, lc};
Point(6+N) = {d, h-a, 0, lc};


Spline(1) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
Line(2) = {16, 2};
Line(3) = {2, 1};
Line(4) = {1, 3};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};

// triple
Physical Point(20) = {3};
Physical Point(21) = {16};

// normal dir
Physical Line(2) = {1};
Physical Line(3) = {4};
Physical Line(4) = {2};

// dir
Physical Line(5) = {3};
Physical Point(5) = {1, 2};


Physical Surface(100) = {6};
