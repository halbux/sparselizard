#include "h1point.h"

using namespace std;


int h1point::count(int order)
{
    if (order <= 0)
        return 0;
    
    return 1;
}

int h1point::count(int order, int dim, int num)
{
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
}



hierarchicalformfunctioncontainer h1point::evalat(int maxorder) 
{    
	element point("point");
    hierarchicalformfunctioncontainer val("h1", point.gettypenumber());

    polynomial formfunc;
    formfunc.set({{{1.0}}});
    val.set(1,0,0,0,0,0,formfunc);

    return val;
}
