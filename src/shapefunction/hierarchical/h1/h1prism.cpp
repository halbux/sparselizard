#include "h1prism.h"

using namespace std;


int h1prism::count(int order)
{
    if (order <= 0)
        return 0;
    
    int countforline = (order+1);
    int countfortriangle = 0.5*(order+1)*(order+2);
    
    return countforline * countfortriangle;
}

int h1prism::count(int order, int dim, int num)
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
            return order-1;
        // Face based form functions:
        case 2:
            // The first two faces are triangular while the last three are quadrangular:
            if (num < 2)
                return 0.5*(order-2)*(order-1);
            else
                return pow(order-1,2);
        // Volume based form functions:
        case 3:
            return 0.5*(order-2)*pow(order-1,2);
    }
}



hierarchicalformfunctioncontainer h1prism::evalat(int maxorder) 
{    
    element prism("prism");
    hierarchicalformfunctioncontainer val("h1", prism.gettypenumber());

    // Get the node list in every edge and face:
    std::vector<int> nodesinedges = prism.getedgesdefinitionsbasedonnodes();                        
    std::vector<int> nodesinfaces = prism.getfacesdefinitionsbasedonnodes();    
    
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
    // phi := (phi+1)/2; 
    phi = 0.5*(phi+1);
    
    vector<polynomial> lambda(7);
    // lambda0 not used
    lambda[1] = 1.0-ki-eta;
    lambda[2] = ki;
    lambda[3] = eta;
    lambda[4] = 1.0-ki-eta;
    lambda[5] = ki;
    lambda[6] = eta;
    
    vector<polynomial> mu(7);
    // mu0 not used
    mu[1] = 1.0-phi;
    mu[2] = 1.0-phi;
    mu[3] = 1.0-phi;
    mu[4] = phi;
    mu[5] = phi;
    mu[6] = phi;
    
    
    ////////// Defining the vertex based form functions (if any):
    
    // Only order 1 has vertex based form functions:
    if (maxorder >= 1)
    {
        // Loop on all nodes:
        for (int node = 0; node < prism.countnodes(); node++)
        {
            polynomial formfunc = lambda[node+1]*mu[node+1];
            val.set(1,0,node,0,0,0,formfunc);
        }
    }
    
    
    ////////// Defining the edge based form functions (if any):

    // Loop on all edges:
    for (int edge = 0; edge < prism.countedges(); edge++)
    {
        // Loop on all possible orientations:
        for (int orientation = 0; orientation < 2; orientation++)
        {
            // Define the nodes e1 and e2 in edge [e1 e2]:
            int e1 = nodesinedges[2*edge+reorderingtoreferenceedgeorientation[orientation][0]] + 1;
            int e2 = nodesinedges[2*edge+reorderingtoreferenceedgeorientation[orientation][1]] + 1;

            // Defining the Legendre polynomials Ls and L for all required orders:
            vector<polynomial> Ls = legendre::Ls(maxorder, lambda[e1]-lambda[e2],lambda[e1]+lambda[e2]);
            vector<polynomial> L = legendre::L(maxorder, 2*mu[e1]-1);

            for (int i = 0; i <= maxorder-2; i++)
            {
                polynomial formfunc;
                if (prism.ishorizontaledge(edge))
                    formfunc = Ls[i+2]*mu[e1];
                else
                    formfunc = L[i+2]*(lambda[e1]+lambda[e2])*0.5; // A factor 0.5 was missing (lambda[e1] = lambda[e2] for horizontal edges)
                
                val.set(i+2,1,edge,orientation,0,0,formfunc);
            }
        }
    }


    ////////// Defining the face based form functions (if any):

    // Loop on all faces:
    for (int face = 0; face < prism.countfaces(); face++)
    {
        if (prism.istriangularface(face))
        {
            // Loop on all possible orientations:
            for (int orientation = 0; orientation < 6; orientation++)
            {            
                // Define the nodes f1, f2 and f3 in face [f1 f2 f3]:
                int f1 = nodesinfaces[3*face+reorderingtoreferencetriangularfaceorientation[orientation][0]] + 1;
                int f2 = nodesinfaces[3*face+reorderingtoreferencetriangularfaceorientation[orientation][1]] + 1;
                int f3 = nodesinfaces[3*face+reorderingtoreferencetriangularfaceorientation[orientation][2]] + 1;

                // Defining the Legendre polynomials Ls and ls for all required orders:
                vector<polynomial> Ls = legendre::Ls(maxorder, lambda[f1]-lambda[f2],lambda[f1]+lambda[f2]);
                vector<polynomial> l = legendre::l(maxorder, 2*lambda[f3]-1);

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

                            polynomial formfunc = Ls[i+2]*lambda[f3]*l[j]*mu[f1];
                            val.set(order,2,face,orientation,ffindex,0,formfunc);

                            ffindex = ffindex + 1;
                        }
                    }
                }
            }
        }
        else
        {
            int offset = 3*prism.counttriangularfaces();
            int quadfaceindex = face-prism.counttriangularfaces();
            
            // Loop on all possible orientations:
            for (int orientation = 0; orientation < 8; orientation++)
            {
                // Define the nodes f1, f2, f3 and f4  in face [f1 f2 f3 f4]:
                int f1 = nodesinfaces[offset+4*quadfaceindex+reorderingtoreferencequadrangularfaceorientation[orientation][0]] + 1;
                int f2 = nodesinfaces[offset+4*quadfaceindex+reorderingtoreferencequadrangularfaceorientation[orientation][1]] + 1;
                int f3 = nodesinfaces[offset+4*quadfaceindex+reorderingtoreferencequadrangularfaceorientation[orientation][2]] + 1;
                int f4 = nodesinfaces[offset+4*quadfaceindex+reorderingtoreferencequadrangularfaceorientation[orientation][3]] + 1;
                
                int f2star;
                // If mu[f1] = mu[f2]:
                if (f1 <= 3 && f2 <= 3 || f1 > 3 && f2 > 3)
                    f2star = f2;
                else
                    f2star = f4;
                
                // Defining the Legendre polynomials for all required orders:
                vector<polynomial> Ls = legendre::Ls(maxorder, lambda[f1]-lambda[f2star],lambda[f1]+lambda[f2star]);
                vector<polynomial> L = legendre::L(maxorder, 2*mu[f1]-1);

                for (int order = 2; order <= maxorder; order++)
                {
                    int ffindex = 0;
                    for (int i = 0; i <= order-2; i++)
                    {
                        for (int j = 0; j <= order-2; j++)
                        {
                            // Skip the part corresponding to a lower order:
                            if (i != order-2 && j != order-2)
                                continue;

                            polynomial formfunc = Ls[i+2]*L[j+2];
                            val.set(order,2,face,orientation,ffindex,0,formfunc);

                            ffindex = ffindex + 1;
                        }
                    }
                }
            }
        }
    }
    
    
    ////////// Defining the volume based form functions (if any):

    // Defining the Legendre polynomials for all required orders:
    vector<polynomial> Ls = legendre::Ls(maxorder, lambda[1]-lambda[2],lambda[1]+lambda[2]);
    vector<polynomial> l = legendre::l(maxorder, 2*lambda[3]-1);
    vector<polynomial> L = legendre::L(maxorder, 2*mu[1]-1);

    for (int order = 3; order <= maxorder; order++)
    {
        int ffindex = 0;
        for (int i = 0; i <= order-3; i++)
        {
            for (int j = 0; j <= order-3; j++)
            {
                for (int k = 0; k <= order-2; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i+j != order-3 && k != order-2)
                        continue;

                    polynomial formfunc =  Ls[i+2]*lambda[3]*l[j]*L[k+2];
                    val.set(order,3,0,0,ffindex,0,formfunc);

                    ffindex = ffindex + 1;
                }
            }
        }
    }
    
    return val;
}
