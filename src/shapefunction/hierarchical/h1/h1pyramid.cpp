#include "h1pyramid.h"

using namespace std;


int h1pyramid::count(int)
{
    std::cout << "Error in 'h1pyramid' object: shape functions not defined for pyramids" << std::endl;
    abort();  
}

int h1pyramid::count(int, int, int)
{
    std::cout << "Error in 'h1pyramid' object: shape functions not defined for pyramids" << std::endl;
    abort();
}



hierarchicalformfunctioncontainer h1pyramid::evalat(int)
{
    std::cout << "Error in 'h1pyramid' object: shape functions not defined for pyramids" << std::endl;
    abort();
}
