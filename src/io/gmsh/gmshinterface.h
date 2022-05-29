// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef GMSHINTERFACE_H
#define GMSHINTERFACE_H

#include <string>
#include "nodes.h"
#include "elements.h"
#include "physicalregions.h"
#include "densemat.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>   
#include <iomanip>
#include "wallclock.h"
#include "physicalregion.h"
#include "element.h"
#include "mystring.h"
#include "polynomial.h"
#include "iodata.h"
#include "lagrangeformfunction.h"

namespace gmshinterface
{
    // Load the mesh from the API to the 'nodes', 'elements' and 'physicalregions' objects.
    void readfromapi(nodes&, elements&, physicalregions&);
    void readwithapi(std::string name, nodes&, elements&, physicalregions&);
    
    // Load the .msh mesh to the 'nodes', 'elements' and 'physicalregions' objects.
    void readfromfile(std::string name, nodes&, elements&, physicalregions&);
    // Write to .msh mesh format:
    void writetofile(std::string name, nodes&, elements&, physicalregions&, disjointregions&, std::vector<int> physicalregionstowrite);
    
    // Write to .pos format:
    void writetofile(std::string name, iodata datatowrite);

    // Write or append the header of a new view in the .pos file:
    void openview(std::string name, std::string viewname, double timetag, bool overwrite);
    // Write a scalar field to the current view in the .pos format:
    void appendtoview(std::string name, int elementtypenumber, densemat coordx, densemat coordy, densemat coordz, densemat compxinterpolated);
    // Write a vector field to the current view in the .pos format:
    void appendtoview(std::string name, int elementtypenumber, densemat coordx, densemat coordy, densemat coordz, densemat compxinterpolated, densemat compyinterpolated, densemat compzinterpolated);
    // Write poly.size() interpolation schemes in the current view in the .pos file:
    void writeinterpolationscheme(std::string name, std::vector<std::vector<polynomial>> poly);
    // Close the current view in the .pos file:
    void closeview(std::string name);
    
    // GMSH comes with its own element type numbering: we translate it to and from ours:
    int convertgmshelementtypenumber(int gmshtypenumber);
    int converttogmshelementtypenumber(int ourtypenumber);
    
    char getelementidentifierinposformat(int ourtypenumber);
};

#endif
