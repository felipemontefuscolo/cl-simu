// Gmsh project created on Mon Oct 25 15:56:52 2010

N = 20;
R = 1;
A = 0.2;
n = 2;   // modo de oscilacao
lc = R/5;

For s In {0:N-1}
	
	t = 1-s/N;
	theta = 2*Pi*t;
	r = ((2*Cos(2*theta)+1)*(R+A) + (1-Cos(2*theta))*Sqrt((2*R-A)*(2*R-A) - 2*A*A))/3;
	x = r * Cos(theta);
	y = r * Sin(theta);
	Point(s+1) = {x,y,0,lc};

EndFor


Spline(1) = {1:N,1};


Line Loop(2) = {1};
Plane Surface(3) = {2};


Physical Line(2) = {1};
Physical Surface(5) = {3};

// tau = R^2 
