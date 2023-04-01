#include "hcurlline.h"

using namespace std;


int hcurlline::count(int order)
{
    if (order < 0)
        return 0;
    
    return order+1;
}

int hcurlline::count(int order, int dim, int num)
{
    // The 'num' input argument is not required here since all nodes, 
    // edges and faces have the same number of form functions. It is
    // however required for prisms and pyramids.
    
    if (order < 0)
        return 0;
    
    switch (dim)
    {
        // Node based form functions:
        case 0:
            return 0;
        // Edge based form functions:
        case 1:
            return order+1;
        // Face based form functions:
        case 2:
            return 0;
        // Volume based form functions:
        case 3:
            return 0;
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer hcurlline::evalat(int maxorder) 
{    
    element line("line");
    hierarchicalformfunctioncontainer val("hcurl", line.gettypenumber());

    // Get the node list in every edge:
    std::vector<int> nodesinedges = line.getedgesdefinitionsbasedonnodes();                        
    
    // Get for every edge orientation the vector reordering the 
    // nodes to bring the edge to its reference orientation 0.
    std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();


    ////////// Define the 'lambda' and 'sigma' polynomials used in Zaglmayr's thesis:
    
    polynomial ki;
    ki.set({{{}},{{{1.0}}}});
    
    // In Zaglmayr's thesis the reference elements are shifted and deformed.
    // Variable change to correspond to our reference element definition:
    // ki := (ki+1)/2;
    ki = 0.5*(ki+1);
    
    vector<polynomial> lambda(3);
    // lambda0 not used
    lambda[1] = 1.0-ki;
    lambda[2] = ki;
    
    vector<polynomial> sigma(3);
    // sigma0 not used
    sigma[1] = 1.0-ki;
    sigma[2] = ki;
    
    
    ////////// Defining the edge based form functions (if any):

    // Loop on all edges:
    for (int edge = 0; edge < line.countedges(); edge++)
    {
        // Loop on all possible orientations:
        for (int orientation = 0; orientation < 2; orientation++)
        {
            // Define the nodes e1 and e2 in edge [e1 e2]:
            int e1 = nodesinedges[2*edge+reorderingtoreferenceedgeorientation[orientation][0]] + 1;
            int e2 = nodesinedges[2*edge+reorderingtoreferenceedgeorientation[orientation][1]] + 1;

            for (int comp = 0; comp < 3; comp++)
            {
                polynomial formfunc = 0.5*(sigma[e1]-sigma[e2]).derivative(comp)*(lambda[e1]+lambda[e2]);
                val.set(0,1,edge,orientation,0,comp,formfunc);
            }

            // Defining the Legendre polynomials L for all required orders:
            vector<polynomial> L = legendre::L(maxorder+1, sigma[e1]-sigma[e2]);

            for (int i = 0; i <= maxorder-1; i++)
            {
                for (int comp = 0; comp < 3; comp++)
                { 
                    polynomial formfunc = (L[i+2]*(lambda[e1]+lambda[e2])).derivative(comp);
                    val.set(i+1,1,edge,orientation,0,comp,formfunc);
                }
            }
        }
    }

    return val;
}

std::vector<bool> hcurlline::isgradienttype(int maxorder)
{
    std::vector<bool> output(count(maxorder,1,0), true);
    // Only the first one is not of grad type:
    output[0] = false;

    return output;
}
