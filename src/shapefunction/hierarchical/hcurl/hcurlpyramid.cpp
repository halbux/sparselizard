#include "hcurlpyramid.h"

using namespace std;


int hcurlpyramid::count(int order)
{
    logs log;
    log.msg() << "Error in 'hcurlpyramid' object: shape functions not defined for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

int hcurlpyramid::count(int order, int dim, int num)
{
    logs log;
    log.msg() << "Error in 'hcurlpyramid' object: shape functions not defined for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer hcurlpyramid::evalat(int maxorder) 
{    
    logs log;
    log.msg() << "Error in 'hcurlpyramid' object: shape functions not defined for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}
