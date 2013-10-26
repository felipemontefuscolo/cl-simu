// Gmsh project created on Mon Oct 25 15:56:52 2010


N = 90;    // numero de pontos (NUMERO IMPAR!!!!)
R = .5;     // raio
A0 = 0.2;   // primeiro modo
A1 = 0.0;  //  segundo modo
A2 = 0.0;  //  terceiro modo  
A3 = 0.0;  //  terceiro modo
A4 = 0.0;
lc = R/1;  // densidade da malha

If((N+1)*.5 != Ceil((N+1)*.5))
  N=N+1;
EndIf

//Point(1) = {L,0,0,lc};
//Point(N) = {0,h,0,lc};
For s In {0:N-1}
	
	t = s/(N-1)/2;
  
	tet = 2*Pi*t  -  Pi/2;

	z = Cos(tet);

	// legendre polinomials
	r = R*(1 + A0*0.5*(3*z*z - 1)   +   A1*0.5*(5*z*z*z - 3*z)   +  A2*0.125*(35*z^4 - 30*z^2 + 3)   +   A3*0.125*(63*z^5 - 70*z^3 + 15*z)) + A4*(231*z^6-315*z^4+105*z^2-5)/16;

	x = r * Cos(tet);
	y = r * Sin(tet);

	Point(s+1) = {x,y,0,lc};

EndFor

Point(N+1) = {0,0,0,lc};

Spline(1) = {1:(N+1)/2};
Spline(2) = {(N+1)/2:N};
Line(3) = {1, N+1};
Line(4) = {N+1, N};
Line(5) = {N+1, (N+1)/2};


Line Loop(6) = {2, -4, 5};
Plane Surface(7) = {6};
Line Loop(8) = {3, 5, -1};
Plane Surface(9) = {8};

Extrude {{0, 1, 0}, {0, 0, 0}, -Pi/2} {
  Surface{7};
}
Extrude {{0, 1, 0}, {0, 0, 0}, -Pi/2} {
  Surface{9};
}



Coherence;
Transfinite Line {11, 2, 15, 25, 1, 3, 13, 5, 4} = 2 Using Progression 1;
Transfinite Surface {16} = {91, 93, 46} Alternated;
Transfinite Surface {32} = {1, 46, 93} Alternated;


Physical Point(1) = {92};
Physical Line(8) = {4, 3};
Physical Surface(6) = {21, 33};
Physical Surface(5) = {7, 9};
Physical Surface(3) = {16, 32};
Physical Volume(39) = {1, 2};

//-dir_tags 1 triple_tags 2 -interf_tags 3 -solid_tags 4,5,6 -feature_tags 7,8,9


Physical Line(6) = {11, 25};
Physical Line(5) = {2, 1};

Physical Point(8) = {91, 1};




Coherence;
