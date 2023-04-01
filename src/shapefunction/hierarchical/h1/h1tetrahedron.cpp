#include "h1tetrahedron.h"

using namespace std;


int h1tetrahedron::count(int order)
{
    if (order <= 0 || targetdim != -1 && targetdim != 3)
        return 0;
    
    return 1.0/6.0*(order+1)*(order+2)*(order+3);
}

int h1tetrahedron::count(int order, int dim, int num)
{
    if (targetdim != -1)
    {
        if (targetdim == 3 && dim == 3)
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
            return 0.5*(order-2)*(order-1);
        // Volume based form functions:
        case 3:
            return 1.0/6.0*(order-3)*(order-2)*(order-1);
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer h1tetrahedron::evalat(int maxorder) 
{    
    std::string type = "h1";
    if (targetdim != -1)
        type = "h1d"+std::to_string(targetdim);
        
    element tetrahedron("tetrahedron");
    hierarchicalformfunctioncontainer val(type, tetrahedron.gettypenumber());
    
    // Get the node list in every edge and face:
    std::vector<int> nodesinedges = tetrahedron.getedgesdefinitionsbasedonnodes();                        
    std::vector<int> nodesinfaces = tetrahedron.getfacesdefinitionsbasedonnodes();    
    
    // Get for every edge and face orientation the vector reordering the 
    // nodes to bring the edge/face to its reference orientation 0.
    std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();
    std::vector<std::vector<int>> reorderingtoreferencetriangularfaceorientation = orientation::getreorderingtoreferencetriangularfaceorientation();
    
    // Store the form function index in a given order:
    std::vector<int> ffindexes(maxorder+1, 0);


    ////////// Define the 'lambda' and 'sigma' polynomials used in Zaglmayr's thesis:
    
    polynomial ki, eta, phi;
    ki.set({{{}},{{{1.0}}}});
    eta.set({{{},{1.0}}});
    phi.set({{{0.0,1.0}}});
    
    vector<polynomial> lambda(5);
    // lambda0 not used
    lambda[1] = 1.0-ki-eta-phi;
    lambda[2] = ki;
    lambda[3] = eta;
    lambda[4] = phi;

    
    ////////// Defining the vertex based form functions (if any):
    
    // Only order 1 has vertex based form functions:
    if (maxorder >= 1)
    {
        // Loop on all nodes:
        for (int node = 0; node < tetrahedron.countnodes(); node++)
        {
            if (targetdim == -1)
                val.set(1,0,node,0,0,0,lambda[node+1]);
            else
            {
                val.set(1,3,0,0,ffindexes[1],0,lambda[node+1]);
                ffindexes[1]++;
            }
        }
    }
    
    
    ////////// Defining the edge based form functions (if any):

    // Loop on all edges:
    for (int edge = 0; edge < tetrahedron.countedges(); edge++)
    {
        // Loop on all possible orientations:
        for (int orientation = 0; orientation < 2; orientation++)
        {
            // Define the nodes e1 and e2 in edge [e1 e2]:
            int e1 = nodesinedges[2*edge+reorderingtoreferenceedgeorientation[orientation][0]] + 1;
            int e2 = nodesinedges[2*edge+reorderingtoreferenceedgeorientation[orientation][1]] + 1;

            // Defining the Legendre polynomials Ls for all required orders:
            vector<polynomial> Ls = legendre::Ls(maxorder, lambda[e1]-lambda[e2],lambda[e1]+lambda[e2]);

            for (int i = 0; i <= maxorder-2; i++)
            {
                polynomial formfunc = Ls[i+2];
                if (targetdim == -1)
                    val.set(i+2,1,edge,orientation,0,0,formfunc);
                else
                {
                    if (orientation == 0)
                    {
                        val.set(i+2,3,0,orientation,ffindexes[i+2],0,formfunc);
                        ffindexes[i+2]++;
                    }
                }
            }
        }
    }


    ////////// Defining the face based form functions (if any):

    // Loop on all faces:
    for (int face = 0; face < tetrahedron.countfaces(); face++)
    {
        // Loop on all possible orientations:
        for (int orientation = 0; orientation < 6; orientation++)
        {            
            // Define the nodes f1, f2 and f3 in face [f1 f2 f3]:
            int f1 = nodesinfaces[3*face+reorderingtoreferencetriangularfaceorientation[orientation][0]] + 1;
            int f2 = nodesinfaces[3*face+reorderingtoreferencetriangularfaceorientation[orientation][1]] + 1;
            int f3 = nodesinfaces[3*face+reorderingtoreferencetriangularfaceorientation[orientation][2]] + 1;

            // Define lambdaF:
            polynomial lambdaF = lambda[f1]+lambda[f2]+lambda[f3];
            
            // Defining the Legendre polynomials Ls and ls for all required orders:
            vector<polynomial> Ls = legendre::Ls(maxorder, lambda[f1]-lambda[f2],lambda[f1]+lambda[f2]);
            vector<polynomial> ls = legendre::ls(maxorder, 2*lambda[f3]-lambdaF,lambdaF);

            for (int order = 3; order <= maxorder; order++)
            {
                int ffindex = 0;
                for (int i = 0; i <= order-3; i++)
                {
                    for (int j = 0; j <= order-3; j++)
                    {
                        // Skip the part corresponding to a lower order:
                        if (i+j != order-3)
                            continue;

                        polynomial formfunc = Ls[i+2]*lambda[f3]*ls[j];
                        if (targetdim == -1)
                            val.set(order,2,face,orientation,ffindex,0,formfunc);
                        else
                        {
                            if (orientation == 0)
                            {
                                val.set(order,3,0,orientation,ffindexes[order],0,formfunc);
                                ffindexes[order]++;
                            }
                        }

                        ffindex = ffindex + 1;
                    }
                }
            }
        }
    }
    
    
    ////////// Defining the volume based form functions (if any):

    // Defining the Legendre polynomials for all required orders:
    vector<polynomial> Ls = legendre::Ls(maxorder, lambda[1]-lambda[2],lambda[1]+lambda[2]);
    vector<polynomial> ls = legendre::ls(maxorder, 2*lambda[3]-(1-lambda[4]),1-lambda[4]);
    vector<polynomial> l = legendre::l(maxorder, 2*lambda[4]-1);

    for (int order = 4; order <= maxorder; order++)
    {
        int ffindex = 0;
        for (int i = 0; i <= order-4; i++)
        {
            for (int j = 0; j <= order-4; j++)
            {
                for (int k = 0; k <= order-4; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i+j+k != order-4)
                        continue;

                    polynomial formfunc = Ls[i+2]*lambda[3]*ls[j]*lambda[4]*l[k];
                    if (targetdim == -1)
                        val.set(order,3,0,0,ffindex,0,formfunc);
                    else
                    {
                        val.set(order,3,0,0,ffindexes[order],0,formfunc);
                        ffindexes[order]++;
                    }

                    ffindex = ffindex + 1;
                }
            }
        }
    }
    
    return val;
}
