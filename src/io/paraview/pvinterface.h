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
#include "mystring.h"

namespace pvinterface
{
    // Write to .vtk format:
    void writetofile(std::string name, iodata datatowrite);
    void writetofile(std::string name, iodata datatowrite, int timestepindex);
    
    // ParaView comes with its own element type numbering:
    int converttoparaviewelementtypenumber(int ourtypenumber);
    // ParaView comes with its own element node ordering:
    std::vector<int> getnodereordering(int ourtypenumber);
};

#endif
