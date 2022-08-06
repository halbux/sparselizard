// This code simulates the Stokes flow of a viscous fluid (water) in a microvalve.
// Reynold's number rho*v*L/mu is very low so that the inertia terms can be neglected.
// The fluid is considered incompressible.
//
// This code can be validated with analytical solutions of a flow between two
// infinite parallel plates with a thin gap. Only the geometry has to be adapted.


#include "sparselizard.h"


using namespace sl;

mesh createmesh(void);

int main(void)
{
    // Define the physical regions that will be used:
    int support = 1, fluid = 2, pillar = 3, membrane = 4, inlet = 7, outlet = 8;
    
    // Create the geometry and the mesh:    
    mesh mymesh = createmesh();
    
    // Write the mesh for display:
    mymesh.write("microvalve.msh");


    // Define the fluid boundary as the intersection between the solid and the fluid regions:
    int solid = selectunion({pillar,support,membrane});
    int fluidboundary = selectintersection({fluid,solid}, 1);
    
    // Dynamic viscosity of water [Pa.s] and density [kg/m3]:
    double mu = 8.9e-4, rho = 1000;
    
    // Field v is the flow velocity. It uses nodal shape functions "h1" with two components.
    // Field p is the relative pressure.
    field v("h1xy"), p("h1");
    
    // Force the flow velocity to 0 at the solid interface:
    v.setconstraint(fluidboundary);
    
    // Set a relative pressure of 100 Pa at the valve inlet and 0 at the outlet:
    p.setconstraint(inlet, 100);
    p.setconstraint(outlet);
    
    // Use an order 1 interpolation for p and 2 for v on the fluid region:
    p.setorder(fluid, 1);
    v.setorder(fluid, 2);
    
    // Define the weak formulation of the Stokes flow problem.
    // The strong form can be found at https://en.wikipedia.org/wiki/Stokes_flow
    formulation viscousflow;

    viscousflow += integral(fluid, predefinedstokes(dof(v), tf(v), dof(p), tf(p), mu, rho, 0, 0) );	
    
    // Generate, solve and save:
    viscousflow.solve();
    
    // Write the p and v fields to file:
    p.write(fluid, "p.vtk", 1);
    v.write(fluid, "v.vtk", 2);
    
    // Output the flowrate for a unit width:
    double flowrate = (-normal(fluid)*v).integrate(inlet, 4);
    std::cout << std::endl << "Flowrate for a unit width: " << flowrate << " m^3/s" << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (flowrate < 1.4845e-7 && flowrate > 1.4844e-7);
}

mesh createmesh(void)
{
    // Give names to the physical region numbers:
    int support = 1, fluid = 2, pillar = 3, membrane = 4, inlet = 7, outlet = 8;

    // Define the x and y geometrical dimensions:
    double y1 = 15e-6, y2 = 25e-6, y3 = 30e-6;
    double x1 = -100e-6, x2 = -50e-6, x3 = -25e-6, x4 = 25e-6, x5 = 50e-6, x6 = 100e-6;
    
    // Define the mesh finesse:
    int n = 10;
    
    shape q1("quadrangle", support, {x1,0,0, x2,0,0, x2,y1,0, x1,y1,0}, {n,n,n,n});
    shape q2("quadrangle", support, {x5,0,0, x6,0,0, x6,y1,0, x5,y1,0}, {n,n,n,n});
    
    shape q3("quadrangle", fluid, {x2,0,0, x3,0,0, x3,y1,0, x2,y1,0}, {n,n,n,n});
    shape q4("quadrangle", fluid, {x4,0,0, x5,0,0, x5,y1,0, x4,y1,0}, {n,n,n,n});
    
    shape q5("quadrangle", pillar, {x3,0,0, x4,0,0, x4,y1,0, x3,y1,0}, {n,n,n,n});
    
    shape q6("quadrangle", support, {x1,y1,0, x2,y1,0, x2,y2,0, x1,y2,0}, {n,n,n,n});
    shape q7("quadrangle", support, {x5,y1,0, x6,y1,0, x6,y2,0, x5,y2,0}, {n,n,n,n});
    
    shape q8("quadrangle", fluid, {x2,y1,0, x3,y1,0, x3,y2,0, x2,y2,0}, {n,n,n,n});
    shape q9("quadrangle", fluid, {x4,y1,0, x5,y1,0, x5,y2,0, x4,y2,0}, {n,n,n,n});
    shape q10("quadrangle", fluid, {x3,y1,0, x4,y1,0, x4,y2,0, x3,y2,0}, {n,n,n,n});
    
    shape q11("quadrangle", membrane, {x1,y2,0, x2,y2,0, x2,y3,0, x1,y3,0}, {n,n,n,n});
    shape q12("quadrangle", membrane, {x5,y2,0, x6,y2,0, x6,y3,0, x5,y3,0}, {n,n,n,n});
    shape q13("quadrangle", membrane, {x2,y2,0, x3,y2,0, x3,y3,0, x2,y3,0}, {n,n,n,n});
    shape q14("quadrangle", membrane, {x4,y2,0, x5,y2,0, x5,y3,0, x4,y3,0}, {n,n,n,n});
    shape q15("quadrangle", membrane, {x3,y2,0, x4,y2,0, x4,y3,0, x3,y3,0}, {n,n,n,n});
    
    // Get the inlet and outlet lines:
    shape l1 = q3.getsons()[0];
    shape l2 = q4.getsons()[0];
    
    l1.setphysicalregion(inlet);
    l2.setphysicalregion(outlet);
    
    // Provide to the mesh all shapes of interest:
    mesh mymesh({q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,l1,l2});
    
    return mymesh;
}

