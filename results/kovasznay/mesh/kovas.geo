// Gmsh project created on Fri May 27 23:13:19 2011
// [-0.5, 1.0] x [-0.5, 0.5]
//
//
//
Point(1) = {-0.5, -0.5, 0, 1.0};
Point(2) = {-0.5, 0.5, 0, 1.0};
Point(3) = {1, 0.5, 0, 1.0};
Point(4) = {1, -0.5, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

Physical Point(2) = {3, 4};
Physical Line(2) = {2, 1, 4};
Physical Line(1) = {3};
Physical Surface(3) = {6};

Transfinite Line {4, 2} = 9 Using Progression 1; //sugestão: 9 | 17 | 33 | 65
Transfinite Line {1, 3} = 7 Using Progression 1; //sugestão: 7 | 13 | 25 | 49
Transfinite Surface {6} Right;

// h = L/(N-1);

// meshes in paper:
// 19x13
// 31x21
// 43x29
// 61x41

