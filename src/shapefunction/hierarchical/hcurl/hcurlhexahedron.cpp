#include "hcurlhexahedron.h"

using namespace std;


int hcurlhexahedron::count(int order)
{
    if (order < 0)
        return 0;
    
    return 3*order*order*(order+1) + 6*2*order*(order+1) + 12*(order+1);
}

int hcurlhexahedron::count(int order, int dim, int num)
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
            return 3*order*order*(order+1);
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer hcurlhexahedron::evalat(int maxorder) 
{    
    element hexahedron("hexahedron");
    hierarchicalformfunctioncontainer val("hcurl", hexahedron.gettypenumber());

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
    
    
    ////////// Defining the volume based form functions (if any):

    // Defining the Legendre polynomials L(2*ki-1), L(2*eta-1) and L(2*phi-1) for all required orders:
    vector<polynomial> Ltwokiminusone = legendre::L(maxorder+1, 2*ki-1);
    vector<polynomial> Ltwoetaminusone = legendre::L(maxorder+1, 2*eta-1);
    vector<polynomial> Ltwophiminusone = legendre::L(maxorder+1, 2*phi-1);

    for (int order = 1; order <= maxorder; order++)
    {
        int ffindex = 0;
        for (int i = 0; i <= order-1; i++)
        {
            for (int j = 0; j <= order-1; j++)
            {
                for (int k = 0; k <= order-1; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i != order-1 && j != order-1 && k != order-1)
                        continue;

                    int ffindexbackup = ffindex;
                    for (int comp = 0; comp < 3; comp++)
                    {
                        ffindex = ffindexbackup;
                                
                        // "Type 1":
                        polynomial phiC1 = (Ltwokiminusone[i+2]*Ltwoetaminusone[j+2]*Ltwophiminusone[k+2]).derivative(comp);
                        val.set(order,3,0,0,ffindex,comp,phiC1);

                        ffindex = ffindex + 1;
                        
                        // "Type 2":
                        polynomial phiC2;
                        if (comp != 1)
                            phiC2 =  phiC1;
                        else
                            phiC2 = -phiC1;

                        val.set(order,3,0,0,ffindex,comp,phiC2);

                        ffindex = ffindex + 1;
                        
                        polynomial phiC2pluspc;
                        if (comp == 0)
                            phiC2pluspc =  phiC1;
                        else
                            phiC2pluspc = -phiC1;

                        val.set(order,3,0,0,ffindex,comp,phiC2pluspc);

                        ffindex = ffindex + 1;
                        
                        // "Type 3":
                        if (i == 0)
                        {
                            polynomial formfunc;
                            if (comp != 0)
                                formfunc.set({{{0.0}}});
                            else
                                formfunc = Ltwoetaminusone[j+2]*Ltwophiminusone[k+2];
                            
                            val.set(order,3,0,0,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;
                        }
                        if (j == 0)
                        {
                            polynomial formfunc;
                            if (comp != 1)
                                formfunc.set({{{0.0}}});
                            else
                                formfunc = Ltwokiminusone[i+2]*Ltwophiminusone[k+2];
                            
                            val.set(order,3,0,0,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;
                        }
                        if (k == 0)
                        {
                            polynomial formfunc;
                            if (comp != 2)
                                formfunc.set({{{0.0}}});
                            else
                                formfunc = Ltwokiminusone[i+2]*Ltwoetaminusone[j+2];
                            
                            val.set(order,3,0,0,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;
                        }
                    }
                }
            }
        }
    }
    
    return val;
}

std::vector<bool> hcurlhexahedron::isgradienttype(int maxorder)
{
    std::vector<bool> output(count(maxorder,3,0), false);

    int ffindex = 0;
    for (int order = 1; order <= maxorder; order++)
    {
        for (int i = 0; i <= order-1; i++)
        {
            for (int j = 0; j <= order-1; j++)
            {
                for (int k = 0; k <= order-1; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i != order-1 && j != order-1 && k != order-1)
                        continue;

                    // "Type 1":
                    output[ffindex] = true; ffindex++;
                    // "Type 2":
                    ffindex++;
                    ffindex++;
                    // "Type 3":
                    if (i == 0)
                        ffindex++;
                    if (j == 0)
                        ffindex++;
                    if (k == 0)
                        ffindex++;
                }
            }
        }
    }

    return output;
}
