// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef NASTRANINTERFACE_H
#define NASTRANINTERFACE_H

#include <string>
#include "nodes.h"
#include "elements.h"
#include "physicalregions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "physicalregion.h"
#include "element.h"
#include "nasdataline.h"

namespace nastraninterface
{
    // Load the mesh to the 'nodes', 'elements' and 'physicalregions' objects.
    void readfromfile(std::string name, nodes&, elements&, physicalregions&);
};

#endif
