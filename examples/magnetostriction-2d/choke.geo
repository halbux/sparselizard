// This mesh file is part of the "Magnetostriction" example of software "GetDP".
// Only minor adaptations were made to include it as a sparselizard example.
// Geometry developed by Mathieu Rossi and Jean Le Besnerais. Copyright (c) 2015 EOMYS



DefineConstant[
  H_sheet = {0.03125, Name "Input/Core parameters/Stack length [m]"},
  nb_sheet = {8, Name "Input/Core parameters/Stack number [integer]"},
  L_sheet = {0.06, Name "Input/Core parameters/Stack width [m]"},
  H_yoke = {L_sheet, Name "Input/Core parameters/Yoke height", ReadOnly 1},
  airgap = {0.00282, Name "Input/Core parameters/Airgap thickness [m]"}
];

DefineConstant[
  square_mesh = {1, Choices{0,1}, Name "Input/Core parameters/Quadrandular mesh"},
  L_yoke = {0.245, Name "Input/Core parameters/Yoke length [m]"},
  winding_thick = {0.03, Name "Input/Winding parameters/winding thickness [m]"}
];

mesh_airgap=6;
mesh_sheet=6;
mesh_column=10;

/////////////////    Yoke    ///////////////////

Point(1) = {-L_yoke/2, 0, 0, 0.004};
Point(2) = {-L_yoke/2+L_sheet, 0, 0, 0.004};
Point(3) = {-L_yoke/2+L_sheet/2, 0, 0, 0.004};
Line(1) = {1, 3};
Line(2) = {3, 2};



Point(4) = {L_yoke/2-L_sheet, 0, 0, 0.004};
Point(5) = {L_yoke/2-L_sheet/2, 0, 0, 0.004};
Point(6) = {L_yoke/2, 0, 0, 0.004};

Line(3) = {4, 5};
Line(4) = {5, 6};

Transfinite Line {1,3} = 15 Using Progression 1.12;
Transfinite Line {2,4} = 15 Using Progression 0.88;

indice={1:4:1};


/////////////// Bouclage pour la colonne /////////////////
For k In {1:nb_sheet}

	If (square_mesh==1)
		nindice[]=Extrude {0, airgap, 0} { Line{indice[]}; Layers{mesh_airgap};Recombine; };
		indice[]={nindice[0],nindice[4],nindice[8],nindice[12]};
		nindice[]=Extrude {0, H_sheet, 0} { Line{indice[]};Layers{mesh_sheet};Recombine;};
		indice[]={nindice[0],nindice[4],nindice[8],nindice[12]};
	EndIf
	If (square_mesh==0)
		nindice[]=Extrude {0, airgap, 0} { Line{indice[]}; };
		indice[]={nindice[0],nindice[4],nindice[8],nindice[12]};
		nindice[]=Extrude {0, H_sheet, 0} { Line{indice[]};};
		indice[]={nindice[0],nindice[4],nindice[8],nindice[12]};
	EndIf


EndFor

If (square_mesh==1)
	nindice[]=Extrude {0, airgap, 0} { Line{indice[]};Layers{mesh_airgap};Recombine; };
EndIf

If (square_mesh==0)
	nindice[]=Extrude {0, airgap, 0} { Line{indice[]}; };
EndIf


indice[]={nindice[0],nindice[4],nindice[8],nindice[12]};
vect_airgap[]={8:8+nb_sheet*32:32,12:12+nb_sheet*32:32,16:16+nb_sheet*32:32,20:20+nb_sheet*32:32};


/////////////// Culasse /////////////////

vect_tole[]={24:24+nb_sheet*32:32,28:28+nb_sheet*32:32,32:32+nb_sheet*32:32,36:36+nb_sheet*32:32};

If (square_mesh==1)
	nindice[]=Extrude {0, H_yoke, 0} { Line{indice[]};Layers{mesh_sheet*2};Recombine;};
	nindice[]=Extrude {0, -H_yoke, 0} { Line{-1,-2,-3,-4};Layers{mesh_sheet*2};Recombine;};
EndIf
If (square_mesh==0)
	nindice[]=Extrude {0, H_yoke, 0} { Line{indice[]};};
	nindice[]=Extrude {0, -H_yoke, 0} { Line{-1,-2,-3,-4};};
EndIf


vect_tole[]={vect_tole[],nindice[1],nindice[6],nindice[11],nindice[16]};
//Printf("Vector tole: %f %f %f %f",) ;

point=(nb_sheet*2+1)*8+1;

Line(500) = {4, 2};
Transfinite Line {500} = 30 Using Progression 1;

If (square_mesh==1)
	nindice[]=Extrude {0, -H_yoke, 0} { Line{500};Layers{mesh_sheet*2};Recombine;};
EndIf
If (square_mesh==0)
	nindice[]=Extrude {0, -H_yoke, 0} { Line{500};};
EndIf


vect_tole[]={vect_tole[],nindice[1]};

Line(550) = {point+1, point+2};
Transfinite Line {550} = 30 Using Progression 1;

If (square_mesh==1)
	nindice[]=Extrude {0, H_yoke, 0} { Line{550};Layers{mesh_sheet*2};Recombine;};
EndIf
If (square_mesh==0)
	nindice[]=Extrude {0, H_yoke, 0} { Line{550};};
EndIf

vect_tole[]={vect_tole[],nindice[1]};



/////////////// conducteur /////////////////

Point(1500) = {-L_yoke/2-0.01, 0.01, 0, 0.01};
Point(1501) = {-L_yoke/2-0.01-winding_thick, 0.01, 0, 0.01};

Point(1502) = {-L_yoke/2+L_sheet+0.01, 0.01, 0, 0.01};
Point(1503) = {-L_yoke/2+L_sheet+0.01+winding_thick, 0.01, 0, 0.01};

Line(1500) = {1500, 1501};
Line(1501) = {1502, 1503};


Point(1504) = {L_yoke/2+0.01, 0.01, 0, 0.01};
Point(1505) = {L_yoke/2+0.01+winding_thick, 0.01, 0, 0.01};


Point(1506) = {L_yoke/2-L_sheet-0.01, 0.01, 0, 0.01};
Point(1507) = {L_yoke/2-L_sheet-0.01-winding_thick, 0.01, 0, 0.01};

Line(1502) = {1504, 1505};
Line(1503) = {1506, 1507};



//Transfinite Line {1500:1503} = 4 Using Progression 1;

nindice[]=Extrude {0, (nb_sheet*H_sheet+(nb_sheet+1)*airgap)-0.02, 0} { Line{-1500,1502};};

//nindice[]=Extrude {0, (nb_sheet*H_sheet+(nb_sheet+1)*airgap)-0.02, 0} { Line{1500,1502};Layers{120};Recombine;};

vect_iplus[]={nindice[1],nindice[6]};
// Why 6 and not 5 as usual ??

nindice[]=Extrude {0, (nb_sheet*H_sheet+(nb_sheet+1)*airgap)-0.02, 0} { Line{1501,-1503};};
//nindice[]=Extrude {0, (nb_sheet*H_sheet+(nb_sheet+1)*airgap)-0.02, 0} { Line{1501,1503};Layers{120};Recombine;};
vect_imoins[]={nindice[1],nindice[5]};


/////////////// Air /////////////////

dist_air=0.1;

h1=H_sheet*nb_sheet+(nb_sheet+1)*airgap+H_yoke+dist_air;
h2=-H_yoke-dist_air;

Point(2516) = {-L_yoke, h1, 0, 1.0};
Point(2517) = {L_yoke, h1, 0, 1.0};
Point(2518) = {L_yoke, h2, 0, 1.0};
Point(2519) = {-L_yoke, h2, 0, 1.0};
Line(2520) = {2516, 2517};
Line(2521) = {2517, 2518};
Line(2522) = {2518, 2519};
Line(2523) = {2519, 2516};

Transfinite Line {2520:2523} = 20 Using Progression 1;
Line Loop(2524) = {2520, 2521, 2522, 2523};


// Ligne interne

vect_interne_gauche[]={11:11+(nb_sheet)*32:16};
vect_culasse_sup_int[]={550};
vect_interne_droitre[]={-14-32*(nb_sheet):-14:16};
vect_culasse_inf_int[]={500};

Line Loop(1550)= {-1503,1518,-1516,-1517};
Line Loop(1551)= {1501,1514,-1512,-1513};



Line Loop(5000) = {vect_interne_gauche[],vect_culasse_sup_int[],vect_interne_droitre[],vect_culasse_inf_int[]};
Plane Surface(6000) = {-5000,1550,1551};

// ligne externe


vect_externe_gauche[]={6:6+(nb_sheet+0.5)*32:16};
vect_culasse_sup_ext[]={5+16*(2*nb_sheet+1),9+16*(2*nb_sheet+1),551,13+16*(2*nb_sheet+1),17+16*(2*nb_sheet+1)};


vect_externe_droitre[]={-19-32*(nb_sheet+0.5):-19:16};
vect_culasse_inf_ext[]={(18+16*(2*nb_sheet+2)),(17+16*(2*nb_sheet+2)),(13+16*(2*nb_sheet+2)),501 ,(9+16*(2*nb_sheet+2)),(5+16*(2*nb_sheet+2)),-(7+16*(2*nb_sheet+2))};
 // 306


Line Loop(6000) = {vect_externe_gauche[],vect_culasse_sup_ext[],vect_externe_droitre[],vect_culasse_inf_ext[]};

Line Loop(1525) = {1504, -1506, 1500, 1505};
Line Loop(1526) = {1508, -1510, -1502, 1509};

Plane Surface(6001) = {-2524,6000,1525,1526};

// Physical

Physical Surface(1)= vect_airgap[];   // airgap
Physical Surface(2)= vect_tole[];       // tole
Physical Surface(3)= vect_iplus[];       // i+
Physical Surface(4)= vect_imoins[];       // i-


Physical Surface(5)= {6000,6001};       // air
Physical Line(6)= {2524};       // air boundaries
Physical Line(7) = {501};   // ground


