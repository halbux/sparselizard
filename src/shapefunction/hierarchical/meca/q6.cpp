#include "q6.h"

using namespace std;


int q6::count(int order)
{
    if (order <= 0)
        return 0;
    
    return 6;
}

int q6::count(int order, int dim, int num)
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
            return 2;
        // Volume based form functions:
        case 3:
            return 0;
    }
}



hierarchicalformfunctioncontainer q6::evalat(int maxorder, vector<double> evaluationpoints) 
{    
	element quadrangle("quadrangle");

    // Reuse the polynomials if available in the universe:
    std::vector<hierarchicalformfunctioncontainer> available = universe::getformfunctionpolys("q6", quadrangle.gettypenumber(), maxorder);
    if (available.size() > 0)
    {
        available[0].evaluate(evaluationpoints);
        return available[0];
    }
    
    hierarchicalformfunctioncontainer val("q6", quadrangle.gettypenumber());

    polynomial ki, eta;
    ki.set({{{}},{{{1.0}}}});
    eta.set({{{},{1.0}}});
    
	// lambda0 not used
	vector<polynomial> lambda(7);
    
    // These form functions are expressed in our reference element and do not need a variable change:
	lambda[5] = 0.5*(1.0-ki*ki);
    lambda[6] = 0.5*(1.0-eta*eta);
    
	// Variable change to correspond to our reference element definition:
	// ki := (ki+1)/2; eta := (eta+1)/2;
    ki = 0.5*(ki+1); eta = 0.5*(eta+1);
    
	lambda[1] = (1.0-ki)*(1.0-eta);
	lambda[2] = ki*(1.0-eta);
    lambda[3] = ki*eta;
    lambda[4] = (1.0-ki)*eta;
    
    for (int i = 0; i < 6; i++)
    {
        if (i < 4)
            val.set(1,0,i,0,0,0,lambda[i+1]);
        else
        {
            for (int orientation = 0; orientation < 8; orientation++)
                val.set(1,2,0,orientation,i-4,0,lambda[i+1]);
        }
    }
    
    universe::setformfunctionpolys("q6", quadrangle.gettypenumber(), maxorder, val);
    
    val.evaluate(evaluationpoints);
    
	return val;
}
