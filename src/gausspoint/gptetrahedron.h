// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef GPTETRAHEDRON_H
#define GPTETRAHEDRON_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gptetrahedron
{   
    int count(int integrationorder); // -1 if not defined
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
