#include "h1pyramid.h"

using namespace std;


int h1pyramid::count(int order)
{
    if (order <= 0)
        return 0;
    
    return 5;// ONLY ORDER 1 IS DEFINED
    
//     int countforline = (order+1);
//     int countfortriangle = 0.5*( (order+1)*(order+1) + (order+1) );
//     
//     return countforline * countfortriangle;
}

int h1pyramid::count(int order, int dim, int num)
{
    
    // The 'num' input argument is required here!
    
    
    if (order <= 0)
        return 0;
    
    switch (dim)
    {
        // Node based form functions:
        case 0:
            return 1;
        // Edge based form functions:
        case 1:
            return 0;//order-1;
        // Face based form functions:
        case 2:
            return 0;
//             // The first two faces are triangular while the last three are quadrangular:
//             if (num < 2)
//                 return 0.5*(order-2)*(order-1);
//             else
//                 return pow(order-1,2);
        // Volume based form functions:
        case 3:
            return 0;
//             return 0.5*(order-2)*pow(order-1,2);
    }
}



hierarchicalformfunctioncontainer h1pyramid::evalat(int maxorder, vector<double> evaluationpoints) 
{    
	element pyramid("pyramid");
    hierarchicalformfunctioncontainer val("h1", pyramid.gettypenumber(), evaluationpoints);

    // Get the node list in every edge and face:
    std::vector<int> nodesinedges = pyramid.getedgesdefinitionsbasedonnodes();						
    std::vector<int> nodesinfaces = pyramid.getfacesdefinitionsbasedonnodes();	
	
	// Get for every edge and face orientation the vector reordering the 
    // nodes to bring the edge/face to its reference orientation 0.
	std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();
	std::vector<std::vector<int>> reorderingtoreferencetriangularfaceorientation = orientation::getreorderingtoreferencetriangularfaceorientation();
	std::vector<std::vector<int>> reorderingtoreferencequadrangularfaceorientation = orientation::getreorderingtoreferencequadrangularfaceorientation();


	////////// Define the 'lambda' and 'mu' polynomials used in Zaglmayr's thesis:
	
    polynomial ki, eta, phi;
    ki.set({{{}},{{{1.0}}}});
    eta.set({{{},{1.0}}});
    phi.set({{{0.0,1.0}}});
    
    // In Zaglmayr's thesis the reference elements are shifted and deformed.
	// Variable change to correspond to our reference element definition:
	// ki := (ki+1)/2; eta := (eta+1)/2;
    ki = 0.5*(ki+1); eta = 0.5*(eta+1);
    
	vector<polynomial> lambda(6);
	// lambda0 not used
	lambda[1] = (1.0-ki)*(1.0-eta)*(1.0-phi);
	lambda[2] = ki*(1.0-eta)*(1.0-phi);
    lambda[3] = ki*eta*(1.0-phi);
    lambda[4] = (1.0-ki)*eta*(1.0-phi);
    lambda[5] = phi;
    

	////////// Defining the vertex based form functions (if any):
	
    // Only order 1 has vertex based form functions:
	if (maxorder >= 1)
	{
		// Loop on all nodes:
		for (int node = 0; node < pyramid.countnodes(); node++)
		{
            val.set(1,0,node,0,0,0,lambda[node+1]);
		}
	}
    
  


    
	return val;
    
}
