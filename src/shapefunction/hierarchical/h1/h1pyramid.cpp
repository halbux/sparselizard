#include "h1pyramid.h"

using namespace std;


int h1pyramid::count(int order)
{
    logs log;
    log.msg() << "Error in 'h1pyramid' object: shape functions not defined for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

int h1pyramid::count(int order, int dim, int num)
{
    logs log;
    log.msg() << "Error in 'h1pyramid' object: shape functions not defined for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer h1pyramid::evalat(int maxorder) 
{    
    logs log;
    log.msg() << "Error in 'h1pyramid' object: shape functions not defined for pyramids" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

