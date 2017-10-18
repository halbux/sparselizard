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
