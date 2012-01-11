// gmsh



lc = 0.5;
R = 1;
x_center = 0;
y_center = 0;
z_center = 0;

Point(1) = {x_center, y_center, z_center, lc};
Point(2) = {x_center - R, y_center, z_center, lc};
Point(4) = {x_center, y_center - R, z_center, lc};
Point(5) = {x_center + R, y_center, z_center, lc};
Point(8) = {x_center, y_center, z_center - R, lc};
Point(11) = {x_center, y_center + R, z_center, lc};
Point(14) = {x_center, y_center, z_center + R, lc};
Circle (1) = {2, 1, 4} Plane{0, 0, 1};
Circle (2) = {4, 1, 5} Plane{0, 0, 1};
Circle (3) = {2, 1, 8} Plane{0, 0, 1};
Circle (4) = {4, 1, 8} Plane{0, 0, 1};
Circle (6) = {2, 1, 11} Plane{0, 0, 1};
Circle (7) = {8, 1, 11} Plane{0, 0, 1};
Circle (9) = {2, 1, 14} Plane{0, 0, 1};
Circle (10) = {11, 1, 14} Plane{0, 0, 1};
Circle (13) = {14, 1, 4} Plane{0, 0, 1};
Circle (15) = {8, 1, 5} Plane{0, 0, 1};
Circle (18) = {11, 1, 5} Plane{0, 0, 1};
Circle (21) = {14, 1, 5} Plane{0, 0, 1};
Line Loop (1000005) = {1, 4, -3};
Ruled Surface (5) = {1000005};
Line Loop (1000008) = {3, 7, -6};
Ruled Surface (8) = {1000008};
Line Loop (1000011) = {6, 10, -9};
Ruled Surface (11) = {1000011};
Line Loop (1000014) = {9, 13, -1};
Ruled Surface (14) = {1000014};
Line Loop (1000017) = {-15, -4, 2};
Ruled Surface (17) = {1000017};
Line Loop (1000020) = {-18, -7, 15};
Ruled Surface (20) = {1000020};
Line Loop (1000023) = {-21, -10, 18};
Ruled Surface (23) = {1000023};
Line Loop (1000026) = {-2, -13, 21};
Ruled Surface (26) = {1000026};

Surface Loop(1000027) = {23, 26, 17, 20, 8, 5, 14, 11};
Volume(1000028) = {1000027};
Physical Surface(2) = {11, 23, 8, 20, 26, 14, 5, 17};
Physical Volume(5) = {1000028};








