// This code simulates in time the laminar, incompressible air flow around a cylinder.
// The forced input air velocity is linearly increased over time from 0 m/s to 0.15 m/s.
// Once the air velocity reaches a threshold value a von Karman vortex street appears.
// The cylinder has a diameter of 14 cm and the truncated air domain is 2 m x 0.8 m.
//
// The Reynolds number (rho*v*D/mu) is 1400 for the max 0.15 m/s velocity.
//
// - rho is the density of air [kg/m3]
// - v is the flow velocity far from the cylinder [m/s]
// - D is the cylinder diameter [m]
// - mu is the dynamic viscosity of air [Pa.s]


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    // Region numbers used in this simulation as defined in the .msh file:
    int fluid = 1, cylinder = 2, inlet = 3, outlet = 4;
    
    // Load the mesh (GMSH format):
    mesh mymesh("channel.msh");
    
    // Define the cylinder skin region:
    int cylskin = regionintersection({fluid,cylinder});
    
    // Field v is the flow velocity. It uses nodal shape functions "h1" with two components in 2D.
    // Field p is the relative pressure.
    field v("h1xy"), p("h1");

    // Force the flow velocity to 0 on the cylinder skin:
    v.setconstraint(cylskin);
    // Force a x-direction flow velocity increasing linearly over time at the inlet:
    v.setconstraint(inlet, array2x1(0.01/6.0*t(),0));
    // Set a 0 relative pressure at the outlet:
    p.setconstraint(outlet);

    // Use an order 1 interpolation for p and 2 for v on the fluid region (satisfies the BB condition):
    p.setorder(fluid, 1); v.setorder(fluid, 2);
    
    // Dynamic viscosity of air [Pa.s] and density [kg/m3] at room temperature and atmospheric pressure:
    double mu = 18e-6, rho = 1.2;
    
    // Define the weak formulation for time-dependent incompressible laminar flow:
    formulation laminarflow;
    
    laminarflow += integral(fluid, predefinednavierstokes(dof(v), tf(v), v, dof(p), tf(p), mu, rho, 0, 0, true) );
    
    // Define the object for an implicit Euler time resolution.
    // The initial field values are taken as the fields are.
    // An all zero initial time derivative is set with 'vec(laminarflow)'.
    impliciteuler eul(laminarflow, vec(laminarflow));
    // Set the relative tolerance on the inner nonlinear iteration:
    eul.settolerance(1e-4);
    
    // Run from 0 sec to 90 sec by steps of 0.2 sec (450 timesteps):
    settime(0);
    for (int i = 0; i < 450; i++)
    {
        // Compute one timestep with an unlimited number of nonlinear iterations (-1):
        eul.next(0.2, -1);
        
        // Write field v with an order 2 interpolation to ParaView .vtu format:
        v.write(fluid, "v" + std::to_string(1000 + i) + ".vtu", 2);
    }
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

