// Gmsh project created on Tue Oct 25 16:06:10 2011

lc = 1;

alpha = 0.5;
beta = 0.7;
w = 1;

k = 1;
a = -5;
b = Acos(alpha/(w*beta))/w + 2*k*Pi;

N=50;
dx = (b-a)/N;

For s In {0:N}
	
	x = a + s*dx;
	y = -alpha*x + (beta-0.00*x*Cos(w*x))*Sin(w*x);
	
	Point(1+s) = {x,y,0,lc};

EndFor

For s In {N+1:2*N-1}
	
	x = a + s*dx;
	z = (2*b-x);
	y = -alpha*z + (beta-0.00*z*Cos(w*z))*Sin(w*z);
	
	Point(1+s) = {x,y,0,lc};

EndFor



Translate {0, -10, 0} {
  Duplicata { Point{1}; }
}
Translate {0, -10, 0} {
  Duplicata { Point{2*N}; }
}

Spline(1) = {1:2*N};

Line(2) = {100, 102};
Line(3) = {102, 101};
Line(4) = {101, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};


Transfinite Line {1, 3} = N+1 Using Progression 1;
Transfinite Line {2, 4} = 2 Using Progression 1;


Transfinite Surface {6};


Recombine Surface {6};
