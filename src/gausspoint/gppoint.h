#ifndef GPPOINT_H
#define GPPOINT_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gppoint
{   
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
