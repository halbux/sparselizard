// sparselizard - Copyright (C) 2017- A. Halbach, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef GPPRISM_H
#define GPPRISM_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gpprism
{   
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
