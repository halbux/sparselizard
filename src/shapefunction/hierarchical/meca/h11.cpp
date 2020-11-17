#include "h11.h"

using namespace std;


int h11::count(int order)
{
    if (order <= 0)
        return 0;
    
    return 11;
}

int h11::count(int order, int dim, [[maybe_unused]] int num)
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
        default:
            return 3;
    }
}



hierarchicalformfunctioncontainer h11::evalat(int maxorder)  // FIXME: use order
{    
    element hexahedron("hexahedron");
    hierarchicalformfunctioncontainer val("h11", hexahedron.gettypenumber());
    
    polynomial ki, eta, phi;
    ki.set({{{}},{{{1.0}}}});
    eta.set({{{},{1.0}}});
    phi.set({{{0.0,1.0}}});

    // lambda0 not used
    vector<polynomial> lambda(12);
    
    // These form functions are expressed in our reference element and do not need a variable change:
    lambda[9] = 0.5*(1.0-ki*ki);
    lambda[10] = 0.5*(1.0-eta*eta);
    lambda[11] = 0.5*(1.0-phi*phi);
    
    // Variable change to correspond to our reference element definition:
    // ki := (ki+1)/2; eta := (eta+1)/2; phi := (phi+1)/2; 
    ki = 0.5*(ki+1); eta = 0.5*(eta+1); phi = 0.5*(phi+1);
    
    lambda[1] = (1.0-ki)*(1.0-eta)*(1.0-phi);
    lambda[2] = ki*(1.0-eta)*(1.0-phi);
    lambda[3] = ki*eta*(1.0-phi);
    lambda[4] = (1.0-ki)*eta*(1.0-phi);
    lambda[5] = (1.0-ki)*(1.0-eta)*phi;
    lambda[6] = ki*(1.0-eta)*phi;
    lambda[7] = ki*eta*phi;
    lambda[8] = (1.0-ki)*eta*phi;
    
    for (int i = 0; i < 11; i++)
    {
        if (i < 8)
            val.set(1,0,i,0,0,0,lambda[i+1]);
        else
            val.set(1,3,0,0,i-8,0,lambda[i+1]);
    }
    
    return val;
}
