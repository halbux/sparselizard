// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// Supported data output formats are:
//
// - GMSH .pos
// - ParaView .vtk


#ifndef IOINTERFACE_H
#define IOINTERFACE_H

#include <string>
#include "iodata.h"

#include "gmshinterface.h"
#include "pvinterface.h"

namespace iointerface
{
    void writetofile(std::string filename, iodata datatowrite, std::string appendtofilename = "");
    // The file format might allow only isoparametric elements:
    bool isonlyisoparametric(std::string filename);
};

#endif
