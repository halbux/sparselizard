#include "h1hexahedron.h"

using namespace std;


int h1hexahedron::count(int order)
{
    if (order <= 0)
        return 0;
    
    return (order+1)*(order+1)*(order+1);
}

int h1hexahedron::count(int order, int dim, int num)
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
            return order-1;
        // Face based form functions:
        case 2:
            return pow(order-1,2);
        // Volume based form functions:
        case 3:
            return pow(order-1,3);
    }
}



hierarchicalformfunctioncontainer h1hexahedron::evalat(int maxorder) 
{    
    element hexahedron("hexahedron");
    hierarchicalformfunctioncontainer val("h1", hexahedron.gettypenumber());

    // Get the node list in every edge and face:
    std::vector<int> nodesinedges = hexahedron.getedgesdefinitionsbasedonnodes();                        
    std::vector<int> nodesinfaces = hexahedron.getfacesdefinitionsbasedonnodes();    
    
    // Get for every edge and face orientation the vector reordering the 
    // nodes to bring the edge/face to its reference orientation 0.
    std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();
    std::vector<std::vector<int>> reorderingtoreferencequadrangularfaceorientation = orientation::getreorderingtoreferencequadrangularfaceorientation();


    ////////// Define the 'lambda' and 'sigma' polynomials used in Zaglmayr's thesis:
    
    polynomial ki, eta, phi;
    ki.set({{{}},{{{1.0}}}});
    eta.set({{{},{1.0}}});
    phi.set({{{0.0,1.0}}});
    
    // In Zaglmayr's thesis the reference elements are shifted and deformed.
    // Variable change to correspond to our reference element definition:
    // ki := (ki+1)/2; eta := (eta+1)/2; phi := (phi+1)/2; 
    ki = 0.5*(ki+1); eta = 0.5*(eta+1); phi = 0.5*(phi+1);
    
    vector<polynomial> lambda(9);
    // lambda0 not used
    lambda[1] = (1.0-ki)*(1.0-eta)*(1.0-phi);
    lambda[2] = ki*(1.0-eta)*(1.0-phi);
    lambda[3] = ki*eta*(1.0-phi);
    lambda[4] = (1.0-ki)*eta*(1.0-phi);
    lambda[5] = (1.0-ki)*(1.0-eta)*phi;
    lambda[6] = ki*(1.0-eta)*phi;
    lambda[7] = ki*eta*phi;
    lambda[8] = (1.0-ki)*eta*phi;
    
    vector<polynomial> sigma(9);
    // sigma0 not used
    sigma[1] = (1.0-ki)+(1.0-eta)+(1.0-phi);
    sigma[2] = ki+(1.0-eta)+(1.0-phi);
    sigma[3] = ki+eta+(1.0-phi);
    sigma[4] = (1.0-ki)+eta+(1.0-phi);
    sigma[5] = (1.0-ki)+(1.0-eta)+phi;
    sigma[6] = ki+(1.0-eta)+phi;
    sigma[7] = ki+eta+phi;
    sigma[8] = (1.0-ki)+eta+phi;
    
    
    ////////// Defining the vertex based form functions (if any):
    
    // Only order 1 has vertex based form functions:
    if (maxorder >= 1)
    {
        // Loop on all nodes:
        for (int node = 0; node < hexahedron.countnodes(); node++)
            val.set(1,0,node,0,0,0,lambda[node+1]);
    }
    
    
    ////////// Defining the edge based form functions (if any):

    // Loop on all edges:
    for (int edge = 0; edge < hexahedron.countedges(); edge++)
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
                val.set(i+2,1,edge,orientation,0,0,formfunc);
            }
        }
    }


    ////////// Defining the face based form functions (if any):

    // Loop on all faces:
    for (int face = 0; face < hexahedron.countfaces(); face++)
    {
        // Loop on all possible orientations:
        for (int orientation = 0; orientation < 8; orientation++)
        {
            // Define the nodes f1, f2, f3 and f4  in face [f1 f2 f3 f4]:
            int f1 = nodesinfaces[4*face+reorderingtoreferencequadrangularfaceorientation[orientation][0]] + 1;
            int f2 = nodesinfaces[4*face+reorderingtoreferencequadrangularfaceorientation[orientation][1]] + 1;
            int f3 = nodesinfaces[4*face+reorderingtoreferencequadrangularfaceorientation[orientation][2]] + 1;
            int f4 = nodesinfaces[4*face+reorderingtoreferencequadrangularfaceorientation[orientation][3]] + 1;

            // Define lambdaF:
            polynomial lambdaF = lambda[f1]+lambda[f2]+lambda[f3]+lambda[f4];
            // Define xiF and etaF:
            polynomial xiF = sigma[f1]-sigma[f2];
            polynomial etaF = sigma[f1]-sigma[f4];

            // Defining the Legendre polynomials L(xiF) and L(etaF) for all required orders:
            vector<polynomial> LxiF = legendre::L(maxorder, xiF);
            vector<polynomial> LetaF = legendre::L(maxorder, etaF);

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

                        polynomial formfunc = LxiF[i+2]*LetaF[j+2]*lambdaF;
                        val.set(order,2,face,orientation,ffindex,0,formfunc);

                        ffindex = ffindex + 1;
                    }
                }
            }
        }
    }
    
    
    ////////// Defining the volume based form functions (if any):

    // Defining the Legendre polynomials L(2*ki-1), L(2*eta-1) and L(2*phi-1) for all required orders:
    vector<polynomial> Ltwokiminusone = legendre::L(maxorder, 2*ki-1);
    vector<polynomial> Ltwoetaminusone = legendre::L(maxorder, 2*eta-1);
    vector<polynomial> Ltwophiminusone = legendre::L(maxorder, 2*phi-1);

    for (int order = 2; order <= maxorder; order++)
    {
        int ffindex = 0;
        for (int i = 0; i <= order-2; i++)
        {
            for (int j = 0; j <= order-2; j++)
            {
                for (int k = 0; k <= order-2; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i != order-2 && j != order-2 && k != order-2)
                        continue;

                    polynomial formfunc = Ltwokiminusone[i+2]*Ltwoetaminusone[j+2]*Ltwophiminusone[k+2];
                    val.set(order,3,0,0,ffindex,0,formfunc);

                    ffindex = ffindex + 1;
                }
            }
        }
    }

    return val;    
}
