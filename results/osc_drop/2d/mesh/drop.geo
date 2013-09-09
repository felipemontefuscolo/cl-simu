// Gmsh project created on Mon Oct 25 15:56:52 2010

N = 40;    // numero de pontos
R = 1;     // raio
A0 = 0.06;   // primeiro modo
A1 = 0.0;  //  segundo modo
A2 = 0.0;  //  terceiro modo  
A3 = 0.0;  //  terceiro modo
A4 = 0.0;
lc = R/30;  // densidade da malha

//Point(1) = {L,0,0,lc};
//Point(N) = {0,h,0,lc};
For s In {0:N-1}
	
	t = s/(N-1)/4;
  
	tet = 2*Pi*t;

	z = Cos(tet);

	// legendre polinomials
	r = R*(1 + A0*0.5*(3*z*z - 1)   +   A1*0.5*(5*z*z*z - 3*z)   +  A2*0.125*(35*z^4 - 30*z^2 + 3)   +   A3*0.125*(63*z^5 - 70*z^3 + 15*z)) + A4*(231*z^6-315*z^4+105*z^2-5)/16;

	x = r * Cos(tet);
	y = r * Sin(tet);

	Point(s+1) = {x,y,0,lc};

EndFor


Point(N+1) = {0,0,0,lc};

Spline(1) = {1:N};

Line(2) = {N, N+1};
Line(3) = {N+1, 1};


Line Loop(4) = {1,2,3};
Plane Surface(5) = -{4};

Rotate {{0, 0, 1}, {0, 0, 0}, 0*Pi/4} {
  Surface{5};
}


Physical Surface(6) = {5};
Physical Line(2) = {1};    // surface
Physical Line(3) = {2};   // solid
Physical Line(4) = {3};   // solid
Physical Point(3) = {N};  // triple
Physical Point(4) = {1};  // triple
Physical Point(5) = {N+1};   // dirichlet



Mesh.Algorithm = 6;



