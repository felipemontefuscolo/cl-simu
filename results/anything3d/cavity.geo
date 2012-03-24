// cavity
cl1 = 0.2;
Point(1) = {0.0, 0.0, 0.0, cl1};
Point(2) = {1.0, 0.0, 0.0, cl1};
Point(3) = {1.0, 1.0, 0.0, cl1};
Point(4) = {0.0, 1.0, 0.0, cl1};

Point(5) = {0.5, 0.0, 0.0, cl1};
Point(6) = {0.0, 0.5, 0.0, cl1};
Point(7) = {0.5, 1.0, 0.0, cl1};
Point(8) = {1.0, 0.5, 0.0, cl1};

Point(9) = {0.5, 0.5, 0.0, cl1};



Extrude {0, 0, 0.5} {
  Point{4, 7, 3, 6, 9, 8, 1, 5, 2};
}
Extrude {0, 0, 0.5} {
  Point{10, 13, 16, 17, 14, 11, 12, 15, 18};
}



Line(19) = {21, 22};
Line(20) = {22, 27};
Line(21) = {27, 26};
Line(22) = {26, 25};
Line(23) = {25, 24};
Line(24) = {24, 19};
Line(25) = {19, 20};
Line(26) = {20, 21};
Line(27) = {20, 23};
Line(28) = {23, 22};
Line(29) = {23, 26};
Line(30) = {23, 24};
Line(31) = {16, 17};
Line(32) = {17, 18};
Line(33) = {18, 15};
Line(34) = {15, 12};
Line(35) = {12, 11};
Line(36) = {11, 10};
Line(37) = {10, 13};
Line(38) = {13, 16};
Line(39) = {13, 14};
Line(40) = {14, 17};
Line(41) = {14, 15};
Line(42) = {14, 11};
Line(43) = {1, 5};
Line(44) = {5, 2};
Line(45) = {2, 8};
Line(46) = {8, 3};
Line(47) = {3, 7};
Line(48) = {7, 4};
Line(49) = {4, 6};
Line(50) = {6, 1};
Line(51) = {6, 9};
Line(52) = {9, 5};
Line(53) = {9, 8};
Line(54) = {9, 7};




Line Loop(55) = {39, 42, 36, 37};
Plane Surface(56) = {55};
Line Loop(57) = {42, -35, -34, -41};
Plane Surface(58) = {57};
Line Loop(59) = {24, 25, 27, 30};
Plane Surface(60) = {-	59};
Line Loop(61) = {23, -30, 29, 22};
Plane Surface(62) = {-61};
Line Loop(63) = {48, 49, 51, 54};
Plane Surface(64) = {-63};
Line Loop(65) = {47, -54, 53, 46};
Plane Surface(66) = {65};
Line Loop(67) = {50, 43, -52, -51};
Plane Surface(68) = {67};
Line Loop(69) = {52, 44, 45, -53};
Plane Surface(70) = {69};
Line Loop(71) = {38, 31, -40, -39};
Plane Surface(72) = {71};
Line Loop(73) = {40, 32, 33, -41};
Plane Surface(74) = {73};
Line Loop(75) = {26, 19, -28, -27};
Plane Surface(76) = {-75};
Line Loop(77) = {29, -21, -20, -28};
Plane Surface(78) = {-77};
Line Loop(79) = {25, -11, -37, 10};
Plane Surface(80) = {-79};
Line Loop(81) = {1, 37, -4, -49};
Plane Surface(82) = {81};
Line Loop(83) = {30, -15, -42, 14};
Plane Surface(84) = {83};
Line Loop(85) = {42, -2, -54, 5};
Plane Surface(86) = {85};
Line Loop(87) = {22, -16, -34, 17};
Plane Surface(88) = {87};
Line Loop(89) = {3, -34, -6, 46};
Plane Surface(90) = {89};
Line Loop(91) = {26, -12, -38, 11};
Plane Surface(92) = {-91};
Line Loop(93) = {38, -7, -50, 4};
Plane Surface(94) = {93};
Line Loop(95) = {28, -13, -40, 14};
Plane Surface(96) = {95};
Line Loop(97) = {40, -8, -52, 5};
Plane Surface(98) = {97};
Line Loop(99) = {21, -17, -33, 18};
Plane Surface(100) = {99};
Line Loop(101) = {33, -6, -45, 9};
Plane Surface(102) = {101};
Line Loop(103) = {19, -13, -31, 12};
Plane Surface(104) = {103};
Line Loop(105) = {31, -8, -43, 7};
Plane Surface(106) = {105};
Line Loop(107) = {27, -14, -39, 11};
Plane Surface(108) = {107};
Line Loop(109) = {39, -5, -51, 4};
Plane Surface(110) = {109};
Line Loop(111) = {24, -10, -36, 15};
Plane Surface(112) = {111};
Line Loop(113) = {36, -1, -48, 2};
Plane Surface(114) = {113};
Line Loop(115) = {20, -18, -32, 13};
Plane Surface(116) = {115};
Line Loop(117) = {32, -9, -44, 8};
Plane Surface(118) = {117};
Line Loop(119) = {29, -17, -41, 14};
Plane Surface(120) = {119};
Line Loop(121) = {41, -6, -53, 5};
Plane Surface(122) = {121};
Line Loop(123) = {23, -15, -35, 16};
Plane Surface(124) = {123};
Line Loop(125) = {35, -2, -47, 3};
Plane Surface(126) = {125};







Surface Loop(127) = {92, 76, 104, 108, 72, 96};
Volume(128) = {127};
Surface Loop(129) = {110, 94, 106, 68, 72, 98};
Volume(130) = {129};
Surface Loop(131) = {108, 80, 60, 112, 84, 56};
Volume(132) = {131};
Surface Loop(133) = {82, 114, 64, 56, 86, 110};
Volume(134) = {133};
Surface Loop(135) = {78, 100, 116, 120, 74, 96};
Volume(136) = {135};
Surface Loop(137) = {98, 118, 102, 70, 74, 122};
Volume(138) = {137};
Surface Loop(139) = {122, 86, 66, 126, 90, 58};
Volume(140) = {139};
Surface Loop(141) = {62, 124, 88, 58, 84, 120};
Volume(142) = {141};




T=2;
Transfinite Line {24, 23, 42, 36, 10, 37, 1, 48, 49, 54, 2, 15, 35, 47, 3, 46, 34, 16, 53, 5, 41, 6, 17, 22, 29, 30, 14, 39, 51, 4, 25, 11, 38, 50, 52, 40, 45, 33, 21, 28, 26, 27, 43, 31, 7, 12, 19, 20, 13, 32, 8, 44, 9, 18} = T Using Progression 1;





Physical Surface(1) = {114, 112, 124, 126};
Physical Surface(2) = {60, 62, 76, 78, 88, 90, 100, 102, 70, 68, 64, 66, 80, 82, 94, 92, 118, 116, 106, 104};



Physical Line(2) = {10, 24, 23, 16, 3, 47, 48, 1};



Physical Volume(9) = {128, 136, 138, 130, 132, 134, 140, 142};



