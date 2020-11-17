#include "hcurlpoint.h"

using namespace std;


int hcurlpoint::count(int)
{
    return 0;
}

int hcurlpoint::count(int, int, int)
{
    return 0;
}



hierarchicalformfunctioncontainer hcurlpoint::evalat(int)
{
    std::cout << "Error in 'hcurlpoint' object: hcurl form functions are not defined for point elements" << std::endl;
    abort();
}
