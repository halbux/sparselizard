#ifndef GPTETRAHEDRON_H
#define GPTETRAHEDRON_H

#include <iostream>
#include <math.h>
#include <vector>
#include "gausspoints.h"

namespace gptetrahedron
{   
    void set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights);
};

#endif
