#include "gppyramid.h"

int gppyramid::count(int integrationorder)
{
    logs log;
    log.msg() << "Error in 'gppyramid' namespace: Gauss points have not been defined yet for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

void gppyramid::set(int integrationorder, std::vector<double>& coordinates, std::vector<double>& weights)
{
    logs log;
    log.msg() << "Error in 'gppyramid' namespace: Gauss points have not been defined yet for pyramids" << std::endl;
    log.error();
}




