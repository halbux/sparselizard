#include "hcurlquadrangle.h"

using namespace std;


int hcurlquadrangle::count(int order)
{
    if (order < 0)
        return 0;
    
    return 2*order*(order+1) + 4*(order+1);
}

int hcurlquadrangle::count(int order, int dim, int num)
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
            return 2*order*(order+1);
        // Volume based form functions:
        case 3:
            return 0;
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer hcurlquadrangle::evalat(int maxorder) 
{    
    element quadrangle("quadrangle");
    hierarchicalformfunctioncontainer val("hcurl", quadrangle.gettypenumber());

    // Get the node list in every edge and face:
    std::vector<int> nodesinedges = quadrangle.getedgesdefinitionsbasedonnodes();                        
    std::vector<int> nodesinfaces = quadrangle.getfacesdefinitionsbasedonnodes();    
    
    // Get for every edge and face orientation the vector reordering the 
    // nodes to bring the edge/face to its reference orientation 0.
    std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();
    std::vector<std::vector<int>> reorderingtoreferencequadrangularfaceorientation = orientation::getreorderingtoreferencequadrangularfaceorientation();


    ////////// Define the 'lambda' and 'sigma' polynomials used in Zaglmayr's thesis:
    
    polynomial ki, eta;
    ki.set({{{}},{{{1.0}}}});
    eta.set({{{},{1.0}}});
    
    // In Zaglmayr's thesis the reference elements are shifted and deformed.
    // Variable change to correspond to our reference element definition:
    // ki := (ki+1)/2; eta := (eta+1)/2;
    ki = 0.5*(ki+1); eta = 0.5*(eta+1);
    
    vector<polynomial> lambda(5);
    // lambda0 not used
    lambda[1] = (1.0-ki)*(1.0-eta);
    lambda[2] = ki*(1.0-eta);
    lambda[3] = ki*eta;
    lambda[4] = (1.0-ki)*eta;
    
    vector<polynomial> sigma(5);
    // sigma0 not used
    sigma[1] = (1.0-ki)+(1.0-eta);
    sigma[2] = ki+(1.0-eta);
    sigma[3] = ki+eta;
    sigma[4] = (1.0-ki)+eta;
    
    
    ////////// Defining the edge based form functions (if any):

    // Loop on all edges:
    for (int edge = 0; edge < quadrangle.countedges(); edge++)
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


    ////////// Defining the face based form functions (if any):

    // Loop on all faces:
    for (int face = 0; face < quadrangle.countfaces(); face++)
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

            // Defining the Legendre polynomials and their derivatives for all required orders:
            vector<polynomial> lxiF = legendre::l(maxorder+1, xiF);
            vector<polynomial> letaF = legendre::l(maxorder+1, etaF);
            vector<polynomial> LxiF = legendre::L(maxorder+1, xiF);
            vector<polynomial> LetaF = legendre::L(maxorder+1, etaF);

            for (int order = 1; order <= maxorder; order++)
            {
                int ffindex = 0;
                for (int i = 0; i <= order-1; i++)
                {
                    for (int j = 0; j <= order-1; j++)
                    {
                        // Skip the part corresponding to a different order:
                        if (i != order-1 && j != order-1)
                            continue;

                        int ffindexbackup = ffindex;
                        for (int comp = 0; comp < 3; comp++)
                        {
                            ffindex = ffindexbackup;
                            
                            // "Type 1":
                            polynomial formfunc = (LxiF[i+2]*LetaF[j+2]*lambdaF).derivative(comp);
                            val.set(order,2,face,orientation,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;

                            // "Type 2":
                            formfunc = ( lxiF[i+2-1]*LetaF[j+2]*xiF.derivative(comp)-LxiF[i+2]*letaF[j+2-1]*etaF.derivative(comp) )*lambdaF;
                            val.set(order,2,face,orientation,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;

                            // "Type 3":
                            if (i == 0)
                            {
                                formfunc = LetaF[j+2]*lambdaF*xiF.derivative(comp);
                                val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                ffindex = ffindex + 1;
                            }
                            if (j == 0)
                            {
                                formfunc = LxiF[i+2]*lambdaF*etaF.derivative(comp);
                                val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                ffindex = ffindex + 1;
                            }
                        }
                    }
                }
            }
        }
    }    

    return val;
}

std::vector<bool> hcurlquadrangle::isgradienttype(int maxorder)
{
    std::vector<bool> output(count(maxorder,2,0), false);

    int ffindex = 0;
    for (int order = 1; order <= maxorder; order++)
    {
        for (int i = 0; i <= order-1; i++)
        {
            for (int j = 0; j <= order-1; j++)
            {
                // Skip the part corresponding to a different order:
                if (i != order-1 && j != order-1)
                    continue;

                // "Type 1":
                output[ffindex] = true; ffindex++;
                // "Type 2":
                ffindex++;
                // "Type 3":
                if (i == 0)
                    ffindex++;
                if (j == 0)
                    ffindex++;
            }
        }
    }

    return output;
}
