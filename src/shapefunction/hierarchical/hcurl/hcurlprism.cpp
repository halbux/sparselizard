#include "hcurlprism.h"

using namespace std;


int hcurlprism::count(int order)
{
    if (order < 0)
        return 0;
    if (order == 0)
        return 9;
    
    return 1.5*(order+1)*order*(order-1) + 2*(order-1)*(order+1) + 3*2*order*(order+1) + 9*(order+1);
}

int hcurlprism::count(int order, int dim, int num)
{
    
    // The 'num' input argument is required here!
    if (order < 0)
        return 0;
    if (order == 0)
    {
        if (dim == 1)
            return 1;
        else
            return 0;
    }
    
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
            // The first two faces are triangular while the last three are quadrangular:
            if (num < 2)
                return (order-1)*(order+1);
            else
                return 2*order*(order+1);
        // Volume based form functions:
        case 3:
            return 1.5*(order+1)*order*(order-1);
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer hcurlprism::evalat(int maxorder) 
{    
    element prism("prism");
    hierarchicalformfunctioncontainer val("hcurl", prism.gettypenumber());

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

            for (int comp = 0; comp < 3; comp++)
            {
                // Lowest order Nedelec function:
                polynomial formfunc;
                if (prism.ishorizontaledge(edge))
                    formfunc = (lambda[e1].derivative(comp)*lambda[e2]-lambda[e1]*lambda[e2].derivative(comp))*mu[e1];
                else
                    formfunc = lambda[e1]*mu[e1].derivative(comp);
                
                val.set(0,1,edge,orientation,0,comp,formfunc);
            }

            // Defining the Legendre polynomials Ls and L for all required orders:
            vector<polynomial> Ls = legendre::Ls(maxorder+1, lambda[e1]-lambda[e2], lambda[e1]+lambda[e2]);
            vector<polynomial> L = legendre::L(maxorder+1, 2*mu[e1]-1);

            for (int i = 0; i <= maxorder-1; i++)
            {
                for (int comp = 0; comp < 3; comp++)
                { 
                    polynomial formfunc;
                    if (prism.ishorizontaledge(edge))
                        formfunc = (Ls[i+2]*mu[e1]).derivative(comp);
                    else
                        formfunc = (L[i+2]*lambda[e1]).derivative(comp);

                    val.set(i+1,1,edge,orientation,0,comp,formfunc);
                }
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

                // Define Ls and l for all required orders:
                vector<polynomial> Ls = legendre::Ls(maxorder, lambda[f1]-lambda[f2], lambda[f1]+lambda[f2]);
                vector<polynomial> l = legendre::l(maxorder, 2*lambda[f3]-1);

                for (int order = 2; order <= maxorder; order++)
                {
                    int ffindex = 0;
                    for (int i = 0; i <= order-2; i++)
                    {
                        for (int j = 0; j <= order-2; j++)
                        {
                            // Skip the part corresponding to a different order:
                            if (i+j != order-2)
                                continue;

                            // Define ui and vj:
                            polynomial ui = Ls[i+2];
                            polynomial vj = lambda[f3]*l[j];

                            int ffindexbackup = ffindex;
                            for (int comp = 0; comp < 3; comp++)
                            {
                                ffindex = ffindexbackup;

                                // "Type 1":
                                polynomial formfunc = (ui*vj*mu[f1]).derivative(comp);
                                val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                ffindex = ffindex + 1;

                                // "Type 2":
                                formfunc = (ui.derivative(comp)*vj-ui*vj.derivative(comp))*mu[f1];
                                val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                ffindex = ffindex + 1;

                                // "Type 3":
                                if (i == 0)
                                {
                                    formfunc = (lambda[f1].derivative(comp)*lambda[f2]-lambda[f1]*lambda[f2].derivative(comp))*vj*mu[f1];
                                    val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                    ffindex = ffindex + 1;
                                }
                            }
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
                int f4 = nodesinfaces[offset+4*quadfaceindex+reorderingtoreferencequadrangularfaceorientation[orientation][3]] + 1;
                
                int f2star;
                // If mu[f1] = mu[f2]:
                if (f1 <= 3 && f2 <= 3 || f1 > 3 && f2 > 3)
                    f2star = f2;
                else
                    f2star = f4;
                
                double alpha;
                if (f2 == f2star)
                    alpha = 1;
                else
                    alpha = -1;

                // Defining the Legendre polynomials and their derivatives for all required orders:
                vector<polynomial> Ls = legendre::Ls(maxorder+1, lambda[f1]-lambda[f2star], lambda[f1]+lambda[f2star]);
                vector<polynomial> L = legendre::L(maxorder+1, 2*mu[f1]-1);

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

                            // Define u and w:
                            polynomial ui = Ls[i+2];
                            polynomial uj = Ls[j+2];
                            polynomial wi = L[i+2];
                            polynomial wj = L[j+2];
                        
                            int ffindexbackup = ffindex;
                            for (int comp = 0; comp < 3; comp++)
                            {
                                ffindex = ffindexbackup;

                                // "Type 1":
                                polynomial formfunc;
                                if (f2 == f2star)
                                    formfunc = (ui*wj).derivative(comp);
                                else
                                    formfunc = (uj*wi).derivative(comp);
                                
                                val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                ffindex = ffindex + 1;

                                // "Type 2":
                                if (f2 == f2star)
                                    formfunc = alpha*(ui.derivative(comp)*wj-ui*wj.derivative(comp));
                                else
                                    formfunc = alpha*(uj.derivative(comp)*wi-uj*wi.derivative(comp));
                                    
                                val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                ffindex = ffindex + 1;

                                // "Type 3":
                                if (i == 0)
                                {
                                    if (f2 == f2star)
                                        formfunc = 2*(lambda[f1].derivative(comp)*lambda[f2star]-lambda[f1]*lambda[f2star].derivative(comp))*wj;
                                    else
                                        formfunc = 2*mu[f1].derivative(comp)*uj;
                                        
                                    val.set(order,2,face,orientation,ffindex,comp,formfunc);

                                    ffindex = ffindex + 1;
                                }
                                if (j == 0)
                                {
                                    if (f2 == f2star)
                                        formfunc = 2*mu[f1].derivative(comp)*ui;
                                    else
                                        formfunc = 2*(lambda[f1].derivative(comp)*lambda[f2star]-lambda[f1]*lambda[f2star].derivative(comp))*wi;
                                        
                                    val.set(order,2,face,orientation,ffindex,comp,formfunc);
                                    
                                    ffindex = ffindex + 1;
                                }
                            }
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
    vector<polynomial> L = legendre::L(maxorder+1, 2*mu[1]-1);
    
    for (int order = 2; order <= maxorder; order++)
    {
        int ffindex = 0;
        for (int i = 0; i <= order-2; i++)
        {
            for (int j = 0; j <= order-2; j++)
            {
                for (int k = 0; k <= order-1; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i+j != order-2 && k != order-1)
                        continue;

                    // Define ui, vj and wk:
                    polynomial ui = Ls[i+2];
                    polynomial vj = lambda[3]*l[j];
                    polynomial wk = L[k+2];
                            
                    int ffindexbackup = ffindex;
                    for (int comp = 0; comp < 3; comp++)
                    {
                        ffindex = ffindexbackup;
                                
                        // "Type 1":
                        polynomial phiC1 = (ui*vj*wk).derivative(comp);
                        val.set(order,3,0,0,ffindex,comp,phiC1);

                        ffindex = ffindex + 1;
                        
                        // "Type 2":
                        polynomial phiC2 = ui.derivative(comp)*vj*wk;
                        val.set(order,3,0,0,ffindex,comp,phiC2);

                        ffindex = ffindex + 1;
                        
                        polynomial phiC2pluspc = ui*vj.derivative(comp)*wk;
                        val.set(order,3,0,0,ffindex,comp,phiC2pluspc);

                        ffindex = ffindex + 1;
                        
                        // "Type 3":
                        if (k == 0)
                        {
                            polynomial phiC3;
                            if (comp != 2)
                                phiC3.set({{{0.0}}});
                            else
                                phiC3 = ui*vj;
                            
                            val.set(order,3,0,0,ffindex,comp,phiC3);

                            ffindex = ffindex + 1;
                        }
                        if (i == 0)
                        {
                            polynomial phiC3 = (lambda[1].derivative(comp)*lambda[2]-lambda[1]*lambda[2].derivative(comp))*vj*wk;
                            val.set(order,3,0,0,ffindex,comp,phiC3);

                            ffindex = ffindex + 1;
                        }
                    }
                }
            }
        }
    }
    
    return val;
}

std::vector<bool> hcurlprism::isgradienttype(int maxorder)
{
    std::vector<bool> output(count(maxorder,3,0), false);

    int ffindex = 0;
    for (int order = 2; order <= maxorder; order++)
    {
        for (int i = 0; i <= order-2; i++)
        {
            for (int j = 0; j <= order-2; j++)
            {
                for (int k = 0; k <= order-1; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i+j != order-2 && k != order-1)
                        continue;

                    // "Type 1":
                    output[ffindex] = true; ffindex++;
                    // "Type 2":
                    ffindex++;
                    ffindex++;
                    // "Type 3":
                    if (k == 0)
                        ffindex++;
                    if (i == 0)
                        ffindex++;
                }
            }
        }
    }

    return output;
}
