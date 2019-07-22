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
    std::cout << "Error in 'hcurlpoint' object: hcurl form functions are not defined for point elements" << std::endl;
    abort();
}
