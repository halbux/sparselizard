#include "hcurlpoint.h"

using namespace std;


int hcurlpoint::count(int order)
{
    return 0;
}

int hcurlpoint::count(int order, int dim, int num)
{
    return 0;
}



hierarchicalformfunctioncontainer hcurlpoint::evalat(int maxorder) 
{    
    logs log;
    log.msg() << "Error in 'hcurlpoint' object: hcurl form functions are not defined for point elements" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}
