#include "hcurlpyramid.h"

using namespace std;


int hcurlpyramid::count(int)
{
    std::cout << "Error in 'hcurlpyramid' object: shape functions not defined for pyramids" << std::endl;
    abort();  
}

int hcurlpyramid::count(int, int, int)
{
    std::cout << "Error in 'hcurlpyramid' object: shape functions not defined for pyramids" << std::endl;
    abort();  
}



hierarchicalformfunctioncontainer hcurlpyramid::evalat(int)
{
    std::cout << "Error in 'hcurlpyramid' object: shape functions not defined for pyramids" << std::endl;
    abort();  
}
