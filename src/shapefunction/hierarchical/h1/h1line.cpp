#include "h1line.h"

using namespace std;


int h1line::count(int order)
{
    if (order <= 0 || targetdim != -1 && targetdim != 1)
        return 0;
    
    return (order+1);
}

int h1line::count(int order, int dim, int num)
{
    if (targetdim != -1)
    {
        if (targetdim == 1 && dim == 1)
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
            return order-1;
        // Face based form functions:
        case 2:
            return 0;
        // Volume based form functions:
        case 3:
            return 0;
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer h1line::evalat(int maxorder) 
{    
    std::string type = "h1";
    if (targetdim != -1)
        type = "h1d"+std::to_string(targetdim);
        
    element line("line");
    hierarchicalformfunctioncontainer val(type, line.gettypenumber());

    // Get the node list in every edge:
    std::vector<int> nodesinedges = line.getedgesdefinitionsbasedonnodes();                        
    
    // Get for every edge orientation the vector reordering the 
    // nodes to bring the edge to its reference orientation 0.
    std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();
    
    // Store the form function index in a given order:
    std::vector<int> ffindexes(maxorder+1, 0);


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
    
    
    ////////// Defining the vertex based form functions (if any):
    
    // Only order 1 has vertex based form functions:
    if (maxorder >= 1)
    {
        // Loop on all nodes:
        for (int node = 0; node < line.countnodes(); node++)
        {
            if (targetdim == -1)
                val.set(1,0,node,0,0,0,lambda[node+1]);
            else
            {
                val.set(1,1,0,0,ffindexes[1],0,lambda[node+1]);
                ffindexes[1]++;
            }
        }
    }
    
    
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

            // Defining the Legendre polynomials L for all required orders:
            vector<polynomial> L = legendre::L(maxorder, sigma[e1]-sigma[e2]);

            for (int i = 0; i <= maxorder-2; i++)
            {
                polynomial formfunc = L[i+2]*(lambda[e1]+lambda[e2]);
                if (targetdim == -1)
                    val.set(i+2,1,edge,orientation,0,0,formfunc);
                else
                {
                    if (orientation == 0)
                    {
                        val.set(i+2,1,0,orientation,ffindexes[i+2],0,formfunc);
                        ffindexes[i+2]++;
                    }
                }
            }
        }
    }
    
    return val;
}
