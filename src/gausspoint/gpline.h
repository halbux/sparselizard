#ifndef GPLINE_H
#define GPLINE_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gpline
{   
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
