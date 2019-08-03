// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// Thanks to R. Harouari for the support of the .vtu output format.


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
#include "universe.h"

namespace pvinterface
{
    void writetovtkfile(std::string name, iodata datatowrite);
    void writetovtufile(std::string name, iodata datatowrite);
    void writetovtkfile(std::string name, iodata datatowrite, int timestepindex);
    void writetovtufile(std::string name, iodata datatowrite, int timestepindex);
    
    // ParaView comes with its own element type numbering:
    int converttoparaviewelementtypenumber(int ourtypenumber);
    // ParaView comes with its own element node ordering:
    std::vector<int> getnodereordering(int ourtypenumber);
};

#endif
