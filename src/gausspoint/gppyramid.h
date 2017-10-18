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
