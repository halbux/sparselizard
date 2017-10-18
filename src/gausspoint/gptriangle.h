#ifndef GPTRIANGLE_H
#define GPTRIANGLE_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gptriangle
{   
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
