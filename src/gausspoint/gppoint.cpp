#include "gppoint.h"

int gppoint::count(int integrationorder)
{
    if (integrationorder >= 0)
        return 1;
    
    return -1;
}

void gppoint::set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights)
{
    coordinates = {0,0,0};
    weights = {1};
}




