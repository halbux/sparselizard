// This mesh file is part of the "Cohomology" example of software "GetDP" (Patrick Dular and Christophe Geuzaine).
// Only minor adaptations were made to include it as a sparselizard example so all the credit for this great 
// geometry and mesh go to the above mentioned persons.



lc = 0.02; // Mesh characteristic length
lc2 = lc;
lc3 = 0;
front3d = 0; // Set to 1 if Frontal 3D mesh algorithm is used
nn = (1./lc)/4.; // Mesh subdivisions per turn, used with Frontal 3D

If(front3d == 1)
  Mesh.Algorithm3D = 4; // Frontal 3D
EndIf
Mesh.Optimize = 1;

DefineConstant
[
  turns = {5, Name "Geometry/Number of coil turns"},
  r = {0.11, Name "Geometry/Coil radius"},
  rc = {0.01, Name "Geometry/Coil wire radius"},
  hc = {0.2, Name "Geometry/Coil height"},
  ht = {0.3, Name "Geometry/Tube height"},
  rt1 = {0.081, Name "Geometry/Tube internal radius"},
  rt2 = {0.092, Name "Geometry/Tube external radius"},
  lb = {1, Name "Geometry/Infinite box width"}
  left = {1, Choices{0,1}, Name "Geometry/Terminals on the left?"}
];

// inductor
p = newp;
Point(p)={0, -r, -hc/2, lc};
Point(p+1)={0, -r+rc, -hc/2, lc};
Point(p+2)={0, -r, -hc/2+rc, lc};
Point(p+3)={0, -r-rc, -hc/2, lc};
Point(p+4)={0, -r, -hc/2-rc,lc};
c = newl;
Circle(c) = {p+1,p,p+2};
Circle(c+1) = {p+2,p,p+3};
Circle(c+2) = {p+3,p,p+4};
Circle(c+3) = {p+4,p,p+1};
ll = newll;
Line Loop(ll) = {c,c+1,c+2,c+3};
s = news;
Plane Surface(s) = {ll};
tmp[] = {s};
vol_coil[] = {};
For j In {1:4*turns+(left?2:0)}
If(front3d == 1)
  tmp[] = Extrude { {0,0,hc/turns/4}, {0,0,1} , {0,0,0} , Pi/2}
                  { Surface {tmp[0]}; Layers {nn / 4}; };
EndIf
If(front3d == 0)
  tmp[] = Extrude { {0,0,hc/turns/4}, {0,0,1} , {0,0,0} , Pi/2}
                  { Surface {tmp[0]}; };
EndIf
  vol_coil[] += tmp[1];
EndFor
If(front3d == 1)
tmp[] = Extrude {(left?-1:1)*lb/2, 0, 0} { Surface{tmp[0]}; Layers{nn}; };
EndIf
If(front3d == 0)
tmp[] = Extrude {(left?-1:1)*lb/2, 0, 0} { Surface{tmp[0]}; };
EndIf
vol_coil[] += tmp[1];
out = tmp[0];
If(front3d == 1)
tmp[] = Extrude {-lb/2, 0, 0} { Surface{s}; Layers{nn}; };
EndIf
If(front3d == 0)
tmp[] = Extrude {-lb/2, 0, 0} { Surface{s}; };
EndIf
vol_coil[] += tmp[1];
in = tmp[0];

// tube
p = newp;
Point(p) = {0, 0, -ht/2, lc2};
Point(p+1) = {rt1, 0, -ht/2, lc2};
Point(p+2) = {0, rt1, -ht/2, lc2};
Point(p+3) = {-rt1, 0, -ht/2, lc2};
Point(p+4) = {0, -rt1, -ht/2, lc2};
Point(p+5) = {rt2, 0, -ht/2, lc2};
Point(p+6) = {0, rt2, -ht/2, lc2};
Point(p+7) = {-rt2, 0, -ht/2, lc2};
Point(p+8) = {0, -rt2, -ht/2, lc2};
c = newc;
Circle(c) = {p+1, p, p+2};
Circle(c+1) = {p+2, p, p+3};
Circle(c+2) = {p+3, p, p+4};
Circle(c+3) = {p+4, p, p+1};
Circle(c+4) = {p+5, p, p+6};
Circle(c+5) = {p+6, p, p+7};
Circle(c+6) = {p+7, p, p+8};
Circle(c+7) = {p+8, p, p+5};
ll = newll;
Line Loop(ll) = {c+4, c+5, c+6, c+7, -c, -(c+1), -(c+2), -(c+3)};
s = news;
Plane Surface(s) = {ll};
If(front3d == 1)
tmp[] = Extrude {0,0,ht}{ Surface{s}; Layers{nn}; };
EndIf
If(front3d == 0)
tmp[] = Extrude {0,0,ht}{ Surface{s}; };
EndIf
vol_tube = tmp[1];

// box
p = newp;
Point(p) = {-lb/2,-lb/2,-lb/2, lc3};
Point(p+1) = {lb/2,-lb/2,-lb/2, lc3};
Point(p+2) = {lb/2,lb/2,-lb/2, lc3};
Point(p+3) = {-lb/2,lb/2,-lb/2, lc3};
Point(p+4) = {-lb/2,-lb/2,lb/2, lc3};
Point(p+5) = {lb/2,-lb/2,lb/2, lc3};
Point(p+6) = {lb/2,lb/2,lb/2, lc3};
Point(p+7) = {-lb/2,lb/2,lb/2, lc3};
l = newl;
Line(l) = {p,p+1};
Line(l+1) = {p+1,p+2};
Line(l+2) = {p+2,p+3};
Line(l+3) = {p+3,p};
Line(l+4) = {p+4,p+5};
Line(l+5) = {p+5,p+6};
Line(l+6) = {p+6,p+7};
Line(l+7) = {p+7,p+4};
Line(l+8) = {p, p+4};
Line(l+9) = {p+1, p+5};
Line(l+10) = {p+2, p+6};
Line(l+11) = {p+3, p+7};
ll = newll;
Line Loop(ll) = Boundary {Surface{in}; };
Line Loop(ll+1) = {l+8, -(l+7), -(l+11), l+3};
Line Loop(ll+2) = Boundary {Surface{out}; };
Line Loop(ll+3) = {l+9, l+5, -(l+10), -(l+1)};
Line Loop(ll+4) = {l,l+1,l+2,l+3};
Line Loop(ll+5) = {l+4,l+5,l+6,l+7};
Line Loop(ll+6) = {l+2, l+11, -(l+6), -(l+10)};
Line Loop(ll+7) = {l, l+9, -(l+4), -(l+8)};
s = news;
tmp[] = {ll+1, ll};
If(left)
  tmp[] += ll+2;
EndIf
Plane Surface(s) = tmp[];
tmp[] = {ll+3};
If(!left)
  tmp[] += ll+2;
EndIf
Plane Surface(s+1) = tmp[];
Plane Surface(s+2) = {ll+4};
Plane Surface(s+3) = {ll+5};
Plane Surface(s+4) = {ll+6};
Plane Surface(s+5) = {ll+7};
sl = newsl;
skin_coil[] = CombinedBoundary{ Volume{vol_coil[]}; };
skin_coil[] -= {in, out};
Surface Loop(sl) = {s:s+5,skin_coil[]};
Surface Loop(sl+1) = CombinedBoundary{ Volume{vol_tube[]}; };
v = newv;
Volume(v) = {sl, sl+1};

COIL = 1;
TUBE = 2;
AIR = 3;
SKIN_COIL = 4;
SKIN_TUBE = 5;
IN = 6;
OUT = 7;
INF = 8;
Physical Volume(COIL) = {vol_coil[]};
Physical Volume(TUBE) = {vol_tube[]};
Physical Volume(AIR) = {v};
Physical Surface(SKIN_COIL) = {skin_coil[]};
Physical Surface(SKIN_TUBE) = CombinedBoundary{ Volume{vol_tube[]}; };
Physical Surface(IN) = in;
Physical Surface(OUT) = out;
Physical Surface(INF) = {s:s+5};

