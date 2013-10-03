// Gmsh project created on Wed Aug 07 2013
// Luzia de Menezes Romanetto - ICMC - USP
// This file is intended to generate a mesh for the simulation of 
// bubbles immersed in a fluid inside a tube

dgh = 133e-6;
dw  = 100e-6;
dgv = 300e-6;

//  V terminal (Stokes) = (2/9)*R^2*g*drho/mu_p  = 3.6296e-06

// capillary waves: dt <= 100 * dx^(3/2)


L=1.00;    				// Global scale factor
lc      = 1.e-5*L;				// Factor for the mesh refinement
lc_wall = 1.e-5*L;				// Factor for the mesh refinement
width  = 600e-6;			// Tube width
height = width; //0.03 / 30;			// Tube height
R = dw/2; 				// Radius of the bubbles
N = 1;					// Number of bubbles TIMES 3
//N = N*3;

// Comment : Here you can give the exact coordinates of the central points 
// of the bubbles
// Point(x) = { coord_x , coord_y, coord_z, lc} , x \in 1..N
Cx = 0.4;
Cx_aux = 0;

Point(1) = {height/2,   width/2,  0, lc};
//Point(2) = {2*dgh+1.5*dw,   dgv,  0, lc};
//Point(3) = {1.5*dgh+dw, 2*dgv,0, lc};

//For t In {1:N/3}
//	//Point(t) = {Cx*L, (19.5-(t-1))*L, 0, lc};
//  Point(3*t-2) = {dgh+R,      (2*t-1)*dgv, 0, lc};
//  Point(3*t-1) = {2*dgh+dw+R, (2*t-1)*dgv, 0, lc};
//  Point(3*t  ) = {1.5*dgh+dw, (2*t  )*dgv, 0, lc};
//	
//	If( Cx == 0.4 )
//		Cx_aux = 0.6;
//	EndIf
//	If( Cx == 0.6 )
//		Cx_aux = 0.4;
//	EndIf	
//	Cx = Cx_aux;	
//EndFor


physicalTagWall_l = 1;
physicalTagWall_r = 2;
physicalTagWall_t = 3;
physicalTagWall_b = 4;
physicalTagSurfaceBubble = 5;
physicalTagOutsideBubble = 100;
physicalTagInsideBubble = 200;
physicalTagFixedPoints = 33;

// ---------------------------------------------------------------------
// ----------------------- Automatic part of the file ------------------
// DO NOT CHANGE AFTER HERE
// -------------- Defining the geometry of bubbles ---------------------

// Bubbles

For t In {1:N}
	// bottom point of the bubble
	Translate {0, -R, 0} { Duplicata{ Point{t}; } }
	
	// top point of the bubble
	Translate {0, R, 0} { Duplicata{ Point{t}; } }

	Circle(2*t-1) = {N+2*t-1, t, N+2*t};
EndFor
For t In {1:N}
	Circle(2*t) = {N+2*t, t, N+2*t-1};
	Line Loop(t) = {2*t-1,2*t};
EndFor

Physical Line(physicalTagSurfaceBubble) = {1:2*N};


// --------------- Defining the geometry of tube -----------------------
// Points
Point(3*N+4) = {0, 0, 0, 2.5*lc_wall};
Point(3*N+3) = {width, 0, 0, 2.5*lc_wall};
Point(3*N+2) = {width, height, 0, 2.5*lc_wall};
Point(3*N+1) = {0, height, 0, 2.5*lc_wall};

// Edges
Line(2*N+1) = {3*N+1, 3*N+2};
Line(2*N+2) = {3*N+2, 3*N+3};
Line(2*N+3) = {3*N+3, 3*N+4};
Line(2*N+4) = {3*N+4, 3*N+1};

Line Loop(N+1) = {2*N+1:2*N+4};

Physical Line(physicalTagWall_l) = {2*N+4};
Physical Line(physicalTagWall_r) = {2*N+2};
Physical Line(physicalTagWall_t) = {2*N+3};
Physical Line(physicalTagWall_b) = {2*N+1};

// ------------- Defining the surface mesh -----------------------------
// Fluid outside the bubbles
Plane Surface(1) = {N+1,1:N};
Physical Surface(physicalTagOutsideBubble) = {1};

// Fluid inside the bubbles
For t In {1:N}
	Plane Surface(t+1) = -{t};
EndFor

Physical Surface(physicalTagInsideBubble) = {2:N+1};
Physical Point(physicalTagFixedPoints) = {3*N+1:3*N+4};

Reverse Surface { 1:N+1 };

