/*************************************************
   DEFINITION OF C-MESH FOR NACA 4-DIGIT AIRFOIL (courtesy of Adrien Crovato)
**************************************************/



//// Physical Groups numbering
IntField = 1;
Airfoil = 2;
Downstream = 3;
Upstream = 4;
Wake = 5;

// Numerical method
maxIter = 15;
relTol = 1.e-5;
absTol = 1.e-5;

//// Freestream
M_infty = 0.7;
alpha = 0*Pi/180;




/* PARAMETERS DEFINED BY USER
*****************************/

// GEOMETRY PARAMETERS
//prop = 8/10;
verticalBorder = 5;
angleLeading = 60*Pi/180;
trailingLength = 5;

// MESHING PARAMETERS
nbElementsLeading = 25;
nbElementsChord = 75;
nbElementsBorder = 30;
nbElementsHorizontalTrailing = 30;

progressionLeading = 1.08;
progressionChord = 1.0;
progressionBorder = 1.15;
progressionHorizontalTrailing = 1.15;



/* GRID DEFINITION
******************/

// GENERATION OF AIRFOIL
N = 201; //number of points defining NACA (must be odd)

Point(1) = {1.000000,0.000000,0};
Point(2) = {0.999753,-0.000035,0};
Point(3) = {0.999013,-0.000141,0};
Point(4) = {0.997781,-0.000317,0};
Point(5) = {0.996057,-0.000562,0};
Point(6) = {0.993844,-0.000876,0};
Point(7) = {0.991144,-0.001258,0};
Point(8) = {0.987958,-0.001707,0};
Point(9) = {0.984292,-0.002222,0};
Point(10) = {0.980147,-0.002801,0};
Point(11) = {0.975528,-0.003443,0};
Point(12) = {0.970440,-0.004147,0};
Point(13) = {0.964888,-0.004909,0};
Point(14) = {0.958877,-0.005729,0};
Point(15) = {0.952414,-0.006603,0};
Point(16) = {0.945503,-0.007531,0};
Point(17) = {0.938153,-0.008510,0};
Point(18) = {0.930371,-0.009537,0};
Point(19) = {0.922164,-0.010610,0};
Point(20) = {0.913540,-0.011726,0};
Point(21) = {0.904508,-0.012883,0};
Point(22) = {0.895078,-0.014079,0};
Point(23) = {0.885257,-0.015310,0};
Point(24) = {0.875056,-0.016574,0};
Point(25) = {0.864484,-0.017868,0};
Point(26) = {0.853553,-0.019189,0};
Point(27) = {0.842274,-0.020535,0};
Point(28) = {0.830656,-0.021904,0};
Point(29) = {0.818712,-0.023291,0};
Point(30) = {0.806454,-0.024694,0};
Point(31) = {0.793893,-0.026111,0};
Point(32) = {0.781042,-0.027539,0};
Point(33) = {0.767913,-0.028974,0};
Point(34) = {0.754521,-0.030414,0};
Point(35) = {0.740877,-0.031856,0};
Point(36) = {0.726995,-0.033296,0};
Point(37) = {0.712890,-0.034733,0};
Point(38) = {0.698574,-0.036163,0};
Point(39) = {0.684062,-0.037582,0};
Point(40) = {0.669369,-0.038988,0};
Point(41) = {0.654508,-0.040378,0};
Point(42) = {0.639496,-0.041747,0};
Point(43) = {0.624345,-0.043094,0};
Point(44) = {0.609072,-0.044414,0};
Point(45) = {0.593691,-0.045705,0};
Point(46) = {0.578217,-0.046962,0};
Point(47) = {0.562667,-0.048182,0};
Point(48) = {0.547054,-0.049362,0};
Point(49) = {0.531395,-0.050499,0};
Point(50) = {0.515705,-0.051587,0};
Point(51) = {0.500000,-0.052625,0};
Point(52) = {0.484295,-0.053608,0};
Point(53) = {0.468605,-0.054534,0};
Point(54) = {0.452946,-0.055397,0};
Point(55) = {0.437333,-0.056195,0};
Point(56) = {0.421783,-0.056924,0};
Point(57) = {0.406309,-0.057581,0};
Point(58) = {0.390928,-0.058163,0};
Point(59) = {0.375655,-0.058666,0};
Point(60) = {0.360504,-0.059087,0};
Point(61) = {0.345492,-0.059424,0};
Point(62) = {0.330631,-0.059674,0};
Point(63) = {0.315938,-0.059834,0};
Point(64) = {0.301426,-0.059902,0};
Point(65) = {0.287110,-0.059876,0};
Point(66) = {0.273005,-0.059754,0};
Point(67) = {0.259123,-0.059535,0};
Point(68) = {0.245479,-0.059217,0};
Point(69) = {0.232087,-0.058799,0};
Point(70) = {0.218958,-0.058280,0};
Point(71) = {0.206107,-0.057661,0};
Point(72) = {0.193546,-0.056940,0};
Point(73) = {0.181288,-0.056119,0};
Point(74) = {0.169344,-0.055197,0};
Point(75) = {0.157726,-0.054176,0};
Point(76) = {0.146447,-0.053056,0};
Point(77) = {0.135516,-0.051839,0};
Point(78) = {0.124944,-0.050527,0};
Point(79) = {0.114743,-0.049121,0};
Point(80) = {0.104922,-0.047624,0};
Point(81) = {0.095492,-0.046037,0};
Point(82) = {0.086460,-0.044364,0};
Point(83) = {0.077836,-0.042608,0};
Point(84) = {0.069629,-0.040770,0};
Point(85) = {0.061847,-0.038854,0};
Point(86) = {0.054497,-0.036863,0};
Point(87) = {0.047586,-0.034800,0};
Point(88) = {0.041123,-0.032668,0};
Point(89) = {0.035112,-0.030471,0};
Point(90) = {0.029560,-0.028212,0};
Point(91) = {0.024472,-0.025893,0};
Point(92) = {0.019853,-0.023517,0};
Point(93) = {0.015708,-0.021088,0};
Point(94) = {0.012042,-0.018607,0};
Point(95) = {0.008856,-0.016078,0};
Point(96) = {0.006156,-0.013503,0};
Point(97) = {0.003943,-0.010884,0};
Point(98) = {0.002219,-0.008223,0};
Point(99) = {0.000987,-0.005521,0};
Point(100) = {0.000247,-0.002779,0};
Point(101) = {0.000000,0.000000,0};
Point(102) = {0.000247,0.002779,0};
Point(103) = {0.000987,0.005521,0};
Point(104) = {0.002219,0.008223,0};
Point(105) = {0.003943,0.010884,0};
Point(106) = {0.006156,0.013503,0};
Point(107) = {0.008856,0.016078,0};
Point(108) = {0.012042,0.018607,0};
Point(109) = {0.015708,0.021088,0};
Point(110) = {0.019853,0.023517,0};
Point(111) = {0.024472,0.025893,0};
Point(112) = {0.029560,0.028212,0};
Point(113) = {0.035112,0.030471,0};
Point(114) = {0.041123,0.032668,0};
Point(115) = {0.047586,0.034800,0};
Point(116) = {0.054497,0.036863,0};
Point(117) = {0.061847,0.038854,0};
Point(118) = {0.069629,0.040770,0};
Point(119) = {0.077836,0.042608,0};
Point(120) = {0.086460,0.044364,0};
Point(121) = {0.095492,0.046037,0};
Point(122) = {0.104922,0.047624,0};
Point(123) = {0.114743,0.049121,0};
Point(124) = {0.124944,0.050527,0};
Point(125) = {0.135516,0.051839,0};
Point(126) = {0.146447,0.053056,0};
Point(127) = {0.157726,0.054176,0};
Point(128) = {0.169344,0.055197,0};
Point(129) = {0.181288,0.056119,0};
Point(130) = {0.193546,0.056940,0};
Point(131) = {0.206107,0.057661,0};
Point(132) = {0.218958,0.058280,0};
Point(133) = {0.232087,0.058799,0};
Point(134) = {0.245479,0.059217,0};
Point(135) = {0.259123,0.059535,0};
Point(136) = {0.273005,0.059754,0};
Point(137) = {0.287110,0.059876,0};
Point(138) = {0.301426,0.059902,0};
Point(139) = {0.315938,0.059834,0};
Point(140) = {0.330631,0.059674,0};
Point(141) = {0.345492,0.059424,0};
Point(142) = {0.360504,0.059087,0};
Point(143) = {0.375655,0.058666,0};
Point(144) = {0.390928,0.058163,0};
Point(145) = {0.406309,0.057581,0};
Point(146) = {0.421783,0.056924,0};
Point(147) = {0.437333,0.056195,0};
Point(148) = {0.452946,0.055397,0};
Point(149) = {0.468605,0.054534,0};
Point(150) = {0.484295,0.053608,0};
Point(151) = {0.500000,0.052625,0};
Point(152) = {0.515705,0.051587,0};
Point(153) = {0.531395,0.050499,0};
Point(154) = {0.547054,0.049362,0};
Point(155) = {0.562667,0.048182,0};
Point(156) = {0.578217,0.046962,0};
Point(157) = {0.593691,0.045705,0};
Point(158) = {0.609072,0.044414,0};
Point(159) = {0.624345,0.043094,0};
Point(160) = {0.639496,0.041747,0};
Point(161) = {0.654508,0.040378,0};
Point(162) = {0.669369,0.038988,0};
Point(163) = {0.684062,0.037582,0};
Point(164) = {0.698574,0.036163,0};
Point(165) = {0.712890,0.034733,0};
Point(166) = {0.726995,0.033296,0};
Point(167) = {0.740877,0.031856,0};
Point(168) = {0.754521,0.030414,0};
Point(169) = {0.767913,0.028974,0};
Point(170) = {0.781042,0.027539,0};
Point(171) = {0.793893,0.026111,0};
Point(172) = {0.806454,0.024694,0};
Point(173) = {0.818712,0.023291,0};
Point(174) = {0.830656,0.021904,0};
Point(175) = {0.842274,0.020535,0};
Point(176) = {0.853553,0.019189,0};
Point(177) = {0.864484,0.017868,0};
Point(178) = {0.875056,0.016574,0};
Point(179) = {0.885257,0.015310,0};
Point(180) = {0.895078,0.014079,0};
Point(181) = {0.904508,0.012883,0};
Point(182) = {0.913540,0.011726,0};
Point(183) = {0.922164,0.010610,0};
Point(184) = {0.930371,0.009537,0};
Point(185) = {0.938153,0.008510,0};
Point(186) = {0.945503,0.007531,0};
Point(187) = {0.952414,0.006603,0};
Point(188) = {0.958877,0.005729,0};
Point(189) = {0.964888,0.004909,0};
Point(190) = {0.970440,0.004147,0};
Point(191) = {0.975528,0.003443,0};
Point(192) = {0.980147,0.002801,0};
Point(193) = {0.984292,0.002222,0};
Point(194) = {0.987958,0.001707,0};
Point(195) = {0.991144,0.001258,0};
Point(196) = {0.993844,0.000876,0};
Point(197) = {0.996057,0.000562,0};
Point(198) = {0.997781,0.000317,0};
Point(199) = {0.999013,0.000141,0};
Point(200) = {0.999753,0.000035,0};

//Points numerotation:

numberLeadingEgde = 101 ;
numberLowerSurface = 81 ;
numberUpperSurface = 121 ;


//Distances:
distanceTrailingLowerPoint = 0.908262 ;
distanceLowerUpperPoint = 0.222905 ;
distanceUpperTrailingPoint = 0.908262 ;
Line(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81};
Line(2) = {81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101};
Line(3) = {101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121};
Line(4) = {121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};


/**********************************
	ELLIPTICAL ARCS
 **********************************/
semiVerticalAxes = verticalBorder;
semiHorizontalAxes = verticalBorder;

Point(10001) = {1, -semiVerticalAxes, 0} ;
Point(10002) = {-semiHorizontalAxes*Cos(angleLeading)+1,
                -semiVerticalAxes*Sin(angleLeading), 0} ;
Point(10003) = {-semiVerticalAxes+1,0,0};
Point(10004) = {-semiHorizontalAxes*Cos(angleLeading)+1,
                semiVerticalAxes*Sin(angleLeading), 0} ;
Point(10005) = {1, semiVerticalAxes, 0} ;


Circle(10001) = {10001, 1, 10002};
Circle(10002) = {10002, 1, 10003};
Circle(10003) = {10003, 1, 10004};
Circle(10004) = {10004, 1, 10005};


Line(10005) = {1, 10001};
Line(10006) = {numberLowerSurface,10002};
Line(10007) = {numberLeadingEgde,10003};
Line(10008) = {numberUpperSurface,10004};
Line(10009) = {1, 10005};


Line Loop(10001) = {1, 10006, -10001, -10005};
Line Loop(10002) = {2, 10007, -10002,-10006};
Line Loop(10003) = {3, 10008, -10003, -10007};
Line Loop(10004) = {4, 10009, -10004, -10008};
Plane Surface(10001) = {10001};
Plane Surface(10002) = {10002};
Plane Surface(10003) = {10003};
Plane Surface(10004) = {10004};

/*******************************
	TRAILING PART
 *******************************/

Point(20001) = {trailingLength+1, verticalBorder, 0} ;
Point(20002) = {trailingLength+1, 0, 0} ;
Point(20003) = {trailingLength+1,-verticalBorder, 0} ;

Line(20001) = {10005, 20001};
Line(20002) = {20002, 20001};
Line(20003) = {20002, 20003};
Line(20004) = {10001, 20003};
Line(20005) = {1, 20002};

Line Loop(20001) = {20005, 20002, -20001, -10009};
Line Loop(20002) = {10005, 20004, -20003, -20005};
Plane Surface(20001) = {20001};
Plane Surface(20002) = {20002};

/*************************
	MESHING
 *************************/

nbNodesBorder = nbElementsBorder + 1;
nbNodesChord = nbElementsChord + 1;
nbNodesLeading = nbElementsLeading + 1;
nbNodesHorizontalTrailing = nbElementsHorizontalTrailing +1;

Delete {Surface{1};} // no grid on airfoil

Transfinite Line {1} = nbNodesChord Using Progression 1/progressionChord;
Transfinite Line {10001} = nbNodesChord Using Progression 1;
Transfinite Line {10005,10006} = nbNodesBorder Using Progression progressionBorder;
Transfinite Surface{10001} = {1, 10001, 10002, numberLowerSurface};
Recombine Surface {10001};

Transfinite Line {10006,10007} = nbNodesBorder Using Progression progressionBorder;
Transfinite Line {2,10002} = nbNodesLeading Using Progression 1/progressionLeading;
Transfinite Surface{10002} = {numberLowerSurface, 10002, 10003, numberLeadingEgde};
Recombine Surface {10002};

Transfinite Line {10007,10008} = nbNodesBorder Using Progression progressionBorder;
Transfinite Line {3,10003} = nbNodesLeading Using Progression progressionLeading;
Transfinite Surface{10003} = {numberLeadingEgde, 10003, 10004, numberUpperSurface};
Recombine Surface {10003};

Transfinite Line {4} = nbNodesChord Using Progression progressionChord;
Transfinite Line {10004} = nbNodesChord Using Progression 1;
Transfinite Line {10008,10009} = nbNodesBorder Using Progression progressionBorder;
Transfinite Surface{10004} = {1, numberUpperSurface, 10004, 10005};
Recombine Surface {10004};

Transfinite Line {10009,20002} = nbNodesBorder Using Progression progressionBorder;
Transfinite Line {20001,20005} = nbNodesHorizontalTrailing Using Progression progressionHorizontalTrailing;
Transfinite Surface{20001} = {1, 20002, 20001, 10005};
Recombine Surface {20001};

Transfinite Line {10005,20003} = nbNodesBorder Using Progression progressionBorder;
Transfinite Line {20004,20005} = nbNodesHorizontalTrailing Using Progression progressionHorizontalTrailing;
Transfinite Surface{20002} = {1, 20002, 20003, 10001};
Recombine Surface {20002};

Physical Line(Upstream) = {10001, 10002, 10003, 10004};
Physical Line(Downstream) = {20001, 20002, 20003, 20004};
Physical Line(Wake) = {2005};
Physical Line(Airfoil) = {1, 2, 3, 4};
Physical Surface(IntField) = {10001, 10002, 10003, 10004, 20001, 20002};
