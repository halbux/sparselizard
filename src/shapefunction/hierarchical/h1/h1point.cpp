#include "h1point.h"

using namespace std;


int h1point::count(int order)
{
    if (order <= 0 || targetdim != -1 && targetdim != 0)
        return 0;
    
    return 1;
}

int h1point::count(int order, int dim, int num)
{
    if (targetdim != -1)
    {
        if (targetdim == 0 && dim == 0)
            return count(order);
        else
            return 0;
    }
    
    // The 'num' input argument is not required here since all nodes, 
    // edges and faces have the same number of form functions. It is
    // however required for prisms and pyramids.
    
    if (order <= 0)
        return 0;
    
    switch (dim)
    {
        // Node based form functions:
        case 0:
            return 1;
        // Edge based form functions:
        case 1:
            return 0;
        // Face based form functions:
        case 2:
            return 0;
        // Volume based form functions:
        case 3:
            return 0;
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer h1point::evalat(int maxorder) 
{    
    std::string type = "h1";
    if (targetdim != -1)
        type = "h1d"+std::to_string(targetdim);

    element point("point");
    hierarchicalformfunctioncontainer val(type, point.gettypenumber());

    polynomial formfunc;
    formfunc.set({{{1.0}}});
    val.set(1,0,0,0,0,0,formfunc);

    return val;
}
