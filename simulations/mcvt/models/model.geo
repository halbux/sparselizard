// Gmsh project created on Mon Nov 09 14:34:02 2020
SetFactory("OpenCASCADE");

Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthFactor = 1;
Mesh.CharacteristicLengthMax = 0.005;

// Construction Parameters
BottomIron = 2/1000;
BottomMagnets = 10/1000;
AirgapBottom = 1/1000;
Iron = 6/1000;
AirgapTop = 1/1000;
TopMagnets = 10/1000;
TopIron = 2/1000;

// Iteration Parameters
BarIntervals = 55;
RotIntervals = 55;

// Main Parameters
Radius = 7.5/100;

// Magnet Parameters
MagnetDiameter = 12/1000;
MagnetRadius = MagnetDiameter/2;

For barsAngle In {0:45:BarIntervals}
  For rotAngle In {0:45:RotIntervals}
    //Deletes the current model (all model entities and their associated meshes).
    Delete Model;
    //Deletes all physical groups.
    Delete Physicals;

    counter = 0;


    // Calculate thickness for these magnets https://www.kjmagnetics.com/thickness.calculator.asp
    // FIRST LAYER
    poles = 2;
    pieces = {24, 16, 12};
    magnetR = MagnetRadius;
    magnetD = {6.2/100, 4.6/100, 3/100};
    layerHeight = BottomMagnets;

    // create surface pieces
    //+
    Disk(1) = {0, 0, 0, Radius, Radius};
    c = 0;
    For j In {0:#pieces[]-1}
      For i In {1:pieces[j]}
        angle = 2 * Pi / pieces[j] * i;
        yangle = magnetD[j] * Cos(angle);
        xangle = magnetD[j] * Sin(angle);
        c = c+1;
        Disk(100 + c) = {xangle, yangle, 0, magnetR, magnetR};    
      EndFor 
    EndFor 

    //+
    BooleanFragments{ Surface{1}; Delete; }{ Surface{101:100+c}; Delete; }
    //+
    Extrude {0, 0, layerHeight} {
      Surface{101:101+c};
    }

    Physical Volume("FirstLayerMagnetsU",2) = {};
    Physical Volume("FirstLayerMagnetsD",3) = {};
    cc = 0;
    For j In {0:#pieces[]-1}
      step = pieces[j]/poles/2; // 2poles / 4 steps
      For i In {0:2*poles-1:2}
        Physical Volume("FirstLayerMagnetsU",2) += {cc+step*i+1:cc+step*(i+1)};
        Physical Volume("FirstLayerMagnetsD",3) += {cc+step*(i+1)+1:cc+step*(i+2)};
      EndFor
      cc = cc + pieces[j];
    EndFor
    //+
    Physical Volume("FirstLayer",101) = {1:1+c};
    counter = 1 + c;
    


    // SECOND LAYER
    zpos = layerHeight + AirgapBottom;
    ironPoles = 8;
    layerHeight = Iron;
    ironS = 2/100;
    ironE = 7/100;

    //+
    Disk(1000) = {0, 0, zpos, Radius, Radius};
    For i In {1:ironPoles}

        angle = (Pi / ironPoles) *  (2*i-1) + barsAngle;
        yangleE = ironE * Cos(angle);
        xangleE = ironE * Sin(angle);
        yangleS = ironS * Cos(angle);
        xangleS = ironS * Sin(angle);
        
        angleN = (Pi / ironPoles) *  (2*i) + barsAngle;
        yangleEN = ironE * Cos(angleN);
        xangleEN = ironE * Sin(angleN);
        yangleSN = ironS * Cos(angleN);
        xangleSN = ironS * Sin(angleN);        

        // Disk(1000 + i) = {xangle, yangle, zpos, ironD * Sin(Pi/ ironPoles/2),ironD * Sin(Pi/ ironPoles/2)};
        Point(1000 + 10*i + 2) = {xangleEN, yangleEN, zpos, 10};
        Point(1000 + 10*i + 1) = {xangleSN, yangleSN, zpos, 10};
        Point(1000 + 10*i + 3) = {xangleE, yangleE, zpos, 10};
        Point(1000 + 10*i + 4) = {xangleS, yangleS, zpos, 10};
        For j In {1:3}
          Line(1000 + 10*i + j) = {1000+10*i+j, 1000+10*i+j+1};
          //Transfinite Curve {1000 + 10*i + j} = 10 Using Progression 1;
        EndFor
        Line(1000 + 10*i + 4) = {1000+10*i+4, 1000+10*i+1};
        //Transfinite Curve {1000 + 10*i + 4} = 10 Using Progression 1;
        Curve Loop(1500 + i) = {1001 + 10*i : 1004+ 10*i};

        Plane Surface(1000+i) = {1500+i};
    EndFor 
    //+
    BooleanFragments{ Surface{1000}; Delete; }{ Surface{1001:1000+ironPoles}; Delete; }
    //+
    Extrude {0, 0, layerHeight} {
      Surface{1001:1001+ironPoles};
    }
    //+
    Physical Volume("IronBars", 4) = {counter+1:counter+ironPoles};
    Physical Volume("SecondLayer",102) = {counter+1:counter+ironPoles+1};
    counter = counter + 1 + ironPoles;


    // Third LAYER
    zpos = zpos + layerHeight + AirgapTop;
    poles = 6;
    pieces = {24, 12, 12};
    magnetR = { MagnetRadius, MagnetRadius, MagnetRadius};
    // magnetD = {6.5/100, 4.6/100, 3/100};
    layerHeight = TopMagnets;

    //+
    Disk(3000) = {0, 0, zpos, Radius, Radius};
    c = 0;
    For j In {0:#pieces[]-1}
      extraAngle = 0;
      If (pieces[j] > 13)
        extraAngle = Pi / pieces[j];
      EndIf
      For i In {1:pieces[j]}
        angle = 2 * Pi / pieces[j] * i + rotAngle + extraAngle;
        yangle = magnetD[j] * Cos(angle);
        xangle = magnetD[j] * Sin(angle);

        c = c+1; 
        Disk(3000 + c) = {xangle, yangle, zpos, magnetR[j], magnetR[j]};   
      EndFor 
    EndFor 
    //+
    BooleanFragments{ Surface{3000}; Delete; }{ Surface{3001:3000+c}; Delete; }
    //+
    Extrude {0, 0, layerHeight} {
      Surface{3001:3001+c};
    }


    //+
    Physical Volume("ThirdLayerMagnetsU",6) = {};
    Physical Volume("ThirdLayerMagnetsD",7) = {};
    cc = counter;
    For j In {0:#pieces[]-1}
      step = pieces[j]/poles/2; // 2poles / 4 steps
      For i In {0:2*poles-1:2}
        Physical Volume("ThirdLayerMagnetsU",6) += {cc+step*i+1:cc+step*(i+1)};
        Physical Volume("ThirdLayerMagnetsD",7) += {cc+step*(i+1)+1:cc+step*(i+2)};
      EndFor
      cc = cc + pieces[j];
    EndFor
    //+
    Physical Volume("ThirdLayer",103) = {counter+1:1+counter+c};
    counter = 1 + c;



    zpos = zpos + layerHeight;
    backIronHeight = TopIron;
    //+
    Cylinder(4000) = {0, 0, -backIronHeight, 0, 0, backIronHeight, Radius, 2*Pi};
    //+
    Cylinder(4001) = {0, 0, zpos, 0, 0, backIronHeight, Radius, 2*Pi};


    //+
    Cylinder(6666) = {0, 0, -backIronHeight, 0, 0, zpos+2*backIronHeight, Radius, 2*Pi};
    //+
    BooleanFragments{
      Volume{6666}; Delete;
    } {
      Physical Volume {101:103,1,2} ; Delete;
    }

    Physical Volume("Air", 10) = {4003, 4004};
    //+
    Recursive Delete {
      Volume{4000,4001};
    }

    Physical Volume("BackIronBottom", 1) = {4005};
    Physical Volume("BackIronTop", 20) = {4002};

    //points = 12;
    //c=0;
    //For d In {2/100:Radius:1/100}
    //For j In {0:points}
    //  For i In {1:points}
    //    angle = 2 * Pi / points * i + rotAngle;
    //    yangle = d * Cos(angle);
    //    xangle = d * Sin(angle);
    //    c = c+1;
    //    Point(99000+c) = {xangle, yangle, 1.6/100, 0.001};
    //  EndFor
    //EndFor
    //EndFor

    //+
    //Transfinite Curve {:} = 10 Using Progression 1;
    //Transfinite Volume {61};
    //Recombine Volume {61};
    // HOW TO Transfinite IRON
    // https://gitlab.onelab.info/gmsh/gmsh/-/issues/358
    // CHECKL ?
    Mesh 3;
    Save Sprintf("parts/Parrot%02gA%02g.msh", rotAngle, barsAngle);

  EndFor
EndFor
