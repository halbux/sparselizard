//+
SetFactory("OpenCASCADE");
Box(1) = {-0.5, -0.6, 0, 1, 1, 1};
//+
Field[1] = MathEval;
//+
Field[1].F = "1";
//+
Field[2] = Restrict;
//+
Field[2].SurfacesList = {4};
//+
Field[1].F = "0.001";
//+
Background Field = 1;
//+
Background Field = 2;
//+
Field[1].F = "1";
//+
Field[1].F = "100";
//+
Field[1].F = "1";
//+
Background Field = 1;
//+
Background Field = 2;
//+
Background Field = -1;
//+
Field[1].F = "0.1";
//+
Field[1].F = "0.0001";
//+
Background Field = 2;
//+
Field[1].F = "1";
//+
Field[1].F = "0.5";
//+
Field[1].F = "0.05";
//+
Field[1].F = "0.005";
//+
Field[1].F = "0.0000005";
//+
Field[1].F = "0.05";
//+
Field[1].F = "0.01";
//+
Field[1].F = "0.1";
//+
Field[1].F = "0.0";
//+
Field[1].F = "0.08";
//+
Field[2].SurfacesList = {};
//+
Field[2].VolumesList = {1};
//+
Field[2].SurfacesList = {};
//+
Field[2].VolumesList = {};
//+
Field[2].VolumesList = {1};
//+
Field[1].F = "0.1";
//+
Field[2].SurfacesList = {};
//+
Delete Field [2];
//+
Field[2] = Restrict;
//+
Field[2].VolumesList = {1};
//+
Field[1].F = "1";
//+
Field[1].F = "0";
//+
Field[1].F = "0.1";
