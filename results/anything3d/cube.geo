// gmsh

lc = 1;
N = 33;  // sugestoes 5 | 9 | 17 | 33

Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {1,1,0,lc};
Point(4) = {0,1,0,lc};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};

Transfinite Line{1:4} = N;
Transfinite Surface{6} = {1,2,3,4};

Extrude {0,0,1} {
  Surface{6}; Layers{N-1};
}



Physical Surface(1) = {27};
Physical Surface(2) = {6, 23, 15, 28, 19};
Physical Line(2) = {4, 10, 8, 2, 18, 9, 14, 3, 1, 22, 11, 13};
Physical Volume(33) = {1};