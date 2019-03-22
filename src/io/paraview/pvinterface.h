// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PVINTERFACE_H
#define PVINTERFACE_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "densematrix.h"
#include "iodata.h"
#include "element.h"

namespace pvinterface
{
	// Write to .vtk format:
	void writetofile(std::string name, iodata datatowrite);

    // ParaView comes with its own element type numbering: we translate it from ours:
    int converttoparaviewelementtypenumber(int ourtypenumber);
};

#endif
