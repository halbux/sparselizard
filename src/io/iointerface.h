// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// Supported data output formats are:
//
// - GMSH .pos
// - ParaView .vtk
// - ParaView .vtu
// - ParaView .pvd

#ifndef IOINTERFACE_H
#define IOINTERFACE_H

#include <string>
#include "iodata.h"
#include "petsc.h"
#include "petscvec.h"

#include "gmshinterface.h"
#include "pvinterface.h"

namespace iointerface
{
    void writetofile(std::string filename, iodata datatowrite, std::string appendtofilename = "");
    // The file format might allow only isoparametric elements:
    bool isonlyisoparametric(std::string filename);
    
    // Group .vtu timestep files in a .pvd file:
    void grouptimesteps(std::string filename, std::vector<std::string> filestogroup, std::vector<double> timevals);
    
    // Write first an int vector then a double vector to ASCII or binary format:
    void write(std::string filename, std::vector<int>& intdata, std::vector<double>& doubledata, bool isbinary);
    
    // Read first an int vector then a double vector from ASCII or binary format:
    void load(std::string filename, std::vector<int>& intdata, std::vector<double>& doubledata, bool isbinary);
};

#endif
