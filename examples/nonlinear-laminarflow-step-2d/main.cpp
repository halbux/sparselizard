// This code simulates the laminar, incompressible water flow past a step. 
//
// A parabolic normal flow velocity is forced at the inlet and 
// a zero pressure is imposed at the outlet.
//
// More information for this standard example can be found in:
//
// "Finite element methods for the incompressible Navier-Stokes equations", A. Segal
//
// The step geometry is as follows:
//
//                                    l2
//                   <--------------------------------->
//                    __________________________________
//                   |       /|\
//            l1     |        |
//      <----------->|        |
//      _____________|    h2 \|/        
//               /|\                              OUTLET
//                |
//      INLET     |
//               \|/ h1
//      ________________________________________________


#include "sparselizardbase.h"


using namespace mathop;

// First four arguments give the geometry dimension, last four arguments give the mesh size:
mesh createmesh(double l1, double h1, double l2, double h2, int nl1, int nl2, int nh1, int nh2);

void sparselizard(void)
{	
    // Region numbers used in this simulation:
    int fluid = 1, inlet = 2, outlet = 3, skin = 4;
    // Height of the inlet [mm]:
    double h1 = 1e-3;
    mesh mymesh = createmesh(2e-3, h1, 12e-3, 1e-3, 30, 150, 20, 50);
    
    // Write the mesh for illustration:
    mymesh.write("pipestep.msh");

    // Define the fluid wall (not including the inlet and outlet):
    int wall = regionexclusion(skin, regionunion({inlet,outlet}));

    // Dynamic viscosity of water [Pa.s] and density [kg/m3]:
    double mu = 8.9e-4, rho = 1000;

    // Field v is the flow velocity. It uses nodal shape functions "h1" with two components.
    // Field p is the pressure. Fields x and y are the space coordinates.
    field v("h1xy"), p("h1"), x("x"), y("y");

    // Force the flow velocity to 0 on the wall:
    v.setconstraint(wall);
    // Set a 0 pressure at the outlet:
    p.setconstraint(outlet);

    // Use an order 1 interpolation for p and 2 for v on the fluid region (satisfies the BB condition):
    p.setorder(fluid, 1);
    v.setorder(fluid, 2);

    // Define the weak formulation for incompressible laminar flow:
    formulation laminarflow;

    laminarflow += integral(fluid, predefinedlaminarflow(dof(v), tf(v), v, dof(p), tf(p), mu, rho) );


    // This loop with the above formulation is a Newton iteration:
    int index = 0; double convergence = 1, velocity = 0.1; 
    while (convergence > 1e-10)
    {
        // Slowly increase the velocity for a high Reynolds (1 m/s flow, 1000 Reynolds still converges). 
        if (velocity < 0.299)
            velocity = velocity + 0.008;

        std::cout << "Flow velocity: " << velocity << " m/s" << std::endl;
        // Force the flow velocity at the inlet (quadratic profile w.r.t. the y axis):
        v.setconstraint(inlet, array2x1(velocity*y*(h1-y)/pow(h1*0.5,2) ,0));

        // Get a measure of the solution for convergence evaluation:
        double measuresol = norm(v).integrate(fluid,2);

        // Generate and solve the laminar flow problem then save to the fields:
        solve(laminarflow);

        convergence = std::abs((norm(v).integrate(fluid,2) - measuresol)/norm(v).integrate(fluid,2));
        std::cout << "Relative solution change: " << convergence << std::endl;

        p.write(fluid, "p.vtk");
        v.write(fluid, "v.vtk", 2);

        index++;
    }

    // Compute the flow velocity norm at position (5,1,0) mm in space:
    double vnorm = norm(v).interpolate(fluid, {5e-3, 1e-3, 0.0})[0];

    // Output the input and output flowrate for a unit width:
    double flowratein = (normal(inlet)*v).integrate(inlet, 4);
    double flowrateout = -(normal(outlet)*v).integrate(outlet, 4);
    std::cout << std::endl << "Flowrate in/out for a unit width: " << flowratein << " / " << flowrateout << " m^3/s" << std::endl;

    // Code validation line. Can be removed.
    std::cout << (vnorm*flowrateout < 2.64489e-05 && vnorm*flowrateout > 2.64485e-05);
}

mesh createmesh(double l1, double h1, double l2, double h2, int nl1, int nl2, int nh1, int nh2)
{
    int fluid = 1, inlet = 2, outlet = 3, skin = 4;

    shape qthinleft("quadrangle", fluid, {0,0,0, l1,0,0, l1,h1,0, 0,h1,0}, {nl1,nh1,nl1,nh1});
    shape qthinright("quadrangle", fluid, {l1,0,0, l1+l2,0,0, l1+l2,h1,0, l1,h1,0}, {nl2,nh1,nl2,nh1});
    shape qthick("quadrangle", fluid, {l1,h1,0, l1+l2,h1,0, l1+l2,h1+h2,0, l1,h1+h2,0}, {nl2,nh2,nl2,nh2});

    shape linlet = qthinleft.getsons()[3];
    linlet.setphysicalregion(inlet);
    shape loutlet = shape("union", outlet, {qthick.getsons()[1],qthinright.getsons()[1]});

    mesh mymesh;
    mymesh.regionskin(skin, fluid);
    mymesh.load({qthinleft,qthinright,qthick,linlet,loutlet});

    mymesh.write("channel.msh");

    return mymesh;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

