// Gmsh project created on Wed Aug 07 2013
// Luzia de Menezes Romanetto - ICMC - USP
// This file is intended to generate a mesh for the simulation of 
// bubbles immersed in a fluid inside a tube

L=1.00;    				// Global scale factor
lc = L/10;				// Factor for the mesh refinement
width  = 1*L;			// Tube width
height = 20*L;			// Tube height
R = 0.2*L; 				// Radius of the bubbles
N = 20;					// Number of bubbles in the tube

// Comment : Here you can give the exact coordinates of the central points 
// of the bubbles
// Point(x) = { coord_x , coord_y, coord_z, lc} , x \in 1..N
Cx = 0.4;
Cx_aux = 0;
For t In {1:N}
	Point(t) = {Cx*L, (19.5-(t-1))*L, 0, lc};	
	
	If( Cx == 0.4 )
		Cx_aux = 0.6;
	EndIf
	If( Cx == 0.6 )
		Cx_aux = 0.4;
	EndIf	
	Cx = Cx_aux;	
EndFor

physicalTagWall = 1;
physicalTagSurfaceBubble = 2;
physicalTagOutsideBubble = 100;
physicalTagInsideBubble = 200;

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
Point(3*N+1) = {0, 0, 0, lc};
Point(3*N+2) = {width, 0, 0, lc};
Point(3*N+3) = {width, height, 0, lc};
Point(3*N+4) = {0, height, 0, lc};

// Edges
Line(2*N+1) = {3*N+1, 3*N+2};
Line(2*N+2) = {3*N+2, 3*N+3};
Line(2*N+3) = {3*N+3, 3*N+4};
Line(2*N+4) = {3*N+4, 3*N+1};

Line Loop(N+1) = {2*N+1:2*N+4};

Physical Line(physicalTagWall) = {2*N+1:2*N+4};

// ------------- Defining the surface mesh -----------------------------
// Fluid outside the bubbles
Plane Surface(1) = {N+1,1:N};
Physical Surface(physicalTagOutsideBubble) = {1};

// Fluid inside the bubbles
For t In {1:N}
	Plane Surface(t+1) = {t};
EndFor

Physical Surface(physicalTagInsideBubble) = {2:N+1};
