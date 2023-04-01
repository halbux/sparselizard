#include "hcurltetrahedron.h"

using namespace std;


int hcurltetrahedron::count(int order)
{
    if (order < 0)
        return 0;
    if (order == 0)
        return 6;
    
    return 0.5*(order-2)*(order-1)*(order+1) + 4*(order-1)*(order+1) + 6*(order+1);
}

int hcurltetrahedron::count(int order, int dim, int num)
{
    // The 'num' input argument is not required here since all nodes, 
    // edges and faces have the same number of form functions. It is
    // however required for prisms and pyramids.
    
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
            return (order-1)*(order+1);
        // Volume based form functions:
        case 3:
            return 0.5*(order-2)*(order-1)*(order+1);
    }
    
    throw std::runtime_error(""); // fix return warning
}



hierarchicalformfunctioncontainer hcurltetrahedron::evalat(int maxorder) 
{    
    element tetrahedron("tetrahedron");
    hierarchicalformfunctioncontainer val("hcurl", tetrahedron.gettypenumber());
    
    // Get the node list in every edge and face:
    std::vector<int> nodesinedges = tetrahedron.getedgesdefinitionsbasedonnodes();                        
    std::vector<int> nodesinfaces = tetrahedron.getfacesdefinitionsbasedonnodes();    
    
    // Get for every edge and face orientation the vector reordering the 
    // nodes to bring the edge/face to its reference orientation 0.
    std::vector<std::vector<int>> reorderingtoreferenceedgeorientation = orientation::getreorderingtoreferenceedgeorientation();
    std::vector<std::vector<int>> reorderingtoreferencetriangularfaceorientation = orientation::getreorderingtoreferencetriangularfaceorientation();


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

            for (int comp = 0; comp < 3; comp++)
            {
                polynomial formfunc = lambda[e1].derivative(comp)*lambda[e2]-lambda[e1]*lambda[e2].derivative(comp);
                val.set(0,1,edge,orientation,0,comp,formfunc);
            }

            // Defining the Legendre polynomials Ls for all required orders:
            vector<polynomial> Ls = legendre::Ls(maxorder+1, lambda[e1]-lambda[e2], lambda[e1]+lambda[e2]);

            for (int i = 0; i <= maxorder-1; i++)
            {
                for (int comp = 0; comp < 3; comp++)
                { 
                    polynomial formfunc = Ls[i+2].derivative(comp);
                    val.set(i+1,1,edge,orientation,0,comp,formfunc);
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
            // Define Ls and ls for all required orders:
            vector<polynomial> Ls = legendre::Ls(maxorder, lambda[f1]-lambda[f2], lambda[f1]+lambda[f2]);
            vector<polynomial> ls = legendre::ls(maxorder, 2*lambda[f3]-lambdaF, lambdaF);

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
                        polynomial vj = lambda[f3]*ls[j];

                        int ffindexbackup = ffindex;
                        for (int comp = 0; comp < 3; comp++)
                        {
                            ffindex = ffindexbackup;
                            
                            // "Type 1":
                            polynomial formfunc = (ui*vj).derivative(comp);
                            val.set(order,2,face,orientation,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;

                            // "Type 2":
                            formfunc = ui.derivative(comp)*vj-ui*vj.derivative(comp);
                            val.set(order,2,face,orientation,ffindex,comp,formfunc);

                            ffindex = ffindex + 1;

                            // "Type 3":
                            if (i == 0)
                            {
                                formfunc = (lambda[f1].derivative(comp)*lambda[f2]-lambda[f1]*lambda[f2].derivative(comp))*vj;
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

    // Defining the Legendre polynomials Ls, ls and l for all required orders:
    vector<polynomial> Ls = legendre::Ls(maxorder, lambda[1]-lambda[2], lambda[1]+lambda[2]);
    vector<polynomial> ls = legendre::ls(maxorder, 2*lambda[3]-(1-lambda[4]), 1-lambda[4]);
    vector<polynomial> l = legendre::l(maxorder, 2*lambda[4]-1);

    for (int order = 3; order <= maxorder; order++)
    {
        int ffindex = 0;
        for (int i = 0; i <= order-3; i++)
        {
            for (int j = 0; j <= order-3; j++)
            {
                for (int k = 0; k <= order-3; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i+j+k != order-3)
                        continue;
                    
                    // Define ui, vj and wk:
                    polynomial ui = Ls[i+2];
                    polynomial vj = lambda[3]*ls[j];
                    polynomial wk = lambda[4]*l[k];

                    int ffindexbackup = ffindex;
                    for (int comp = 0; comp < 3; comp++)
                    {
                        ffindex = ffindexbackup;
                                
                        // "Type 1":
                        polynomial phiC1 = (ui*vj*wk).derivative(comp);
                        val.set(order,3,0,0,ffindex,comp,phiC1);

                        ffindex = ffindex + 1;
                        
                        // "Type 2":
                        polynomial phiC2 = ui.derivative(comp)*vj*wk-ui*vj.derivative(comp)*wk+ui*vj*wk.derivative(comp);
                        val.set(order,3,0,0,ffindex,comp,phiC2);

                        ffindex = ffindex + 1;
                        
                        polynomial phiC2pluspc = ui.derivative(comp)*vj*wk+ui*vj.derivative(comp)*wk-ui*vj*wk.derivative(comp);
                        val.set(order,3,0,0,ffindex,comp,phiC2pluspc);

                        ffindex = ffindex + 1;
                        
                        // "Type 3":
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

std::vector<bool> hcurltetrahedron::isgradienttype(int maxorder)
{
    std::vector<bool> output(count(maxorder,3,0), false);

    int ffindex = 0;
    for (int order = 3; order <= maxorder; order++)
    {
        for (int i = 0; i <= order-3; i++)
        {
            for (int j = 0; j <= order-3; j++)
            {
                for (int k = 0; k <= order-3; k++)
                {
                    // Skip the part corresponding to a lower order:
                    if (i+j+k != order-3)
                        continue;
                    
                    // "Type 1":
                    output[ffindex] = true; ffindex++;
                    // "Type 2":
                    ffindex++;
                    ffindex++;
                    // "Type 3":
                    if (i == 0)
                        ffindex++;
                }
            }
        }
    }

    return output;
}
