// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef GPPYRAMID_H
#define GPPYRAMID_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gppyramid
{   
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
