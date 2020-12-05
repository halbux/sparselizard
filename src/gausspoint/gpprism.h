// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef GPPRISM_H
#define GPPRISM_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gpprism
{   
    int count(int integrationorder); // -1 if not defined
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
