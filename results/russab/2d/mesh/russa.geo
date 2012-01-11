// Gmsh project created on Tue Oct 25 16:12:26 2011// Gmsh project created on Tue Oct 25 16:06:10 2011

lc = 0.08;

alpha = 0.5;
beta = 0.7;
w = 1;

k = 2;
a = -5;
b = Acos(alpha/(w*beta))/w + 2*k*Pi;

N=300;
dx = (b-a)/N;

For s In {0:N/10}
	
	x = a + s*dx;
	y = -alpha*x + (beta-0.00*x*Cos(w*x))*Sin(w*x);
	
	Point(1+s) = {x,y,0,lc};

EndFor

Spline(1) = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};





Point(32) = {-3.2, 3.5, 0,   lc};




Ellipse(2) = {1, 17, 17, 32};
Ellipse(3) = {32, 18, 18, 31};
Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};


Physical Surface(100) = {5};
Physical Line(2) = {2, 3};

Physical Point(20) = {1};
Physical Point(21) = {31};


Physical Line(3) = {1};
// Progression
Transfinite Line {1} = 40 Using Bump 0.3;



