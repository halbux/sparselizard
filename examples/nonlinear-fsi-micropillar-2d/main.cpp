// STATUS: almost validated, up to > 30 degrees pillar tilt, match better than 1% on peak deflection
// NEXT STEPS: cleaning code a little + validating with 4 significant digits

// This code simulates the fluid-structure interaction between an incompressible water 
// flow in a microchannel and a polyimide micropillar. The micropillar is modeled as 
// an elastic structure. Small-strain geometric-nonlinearity is taken into account.
// A monolithic fluid-structure coupling is used. 
//
// A parabolic normal flow velocity is forced at the inlet and 
// a zero pressure is imposed at the outlet.
//
//
// The microchannel geometry is as follows (l = 350 um, h = 120 um):
//
//                             l
//      <---------------------------------------------->
//      ________________________________    ____________
//                                      |  |
//               /|\                    |  |      OUTLET
//                |                     |  |
//                |                     |  |
//                |       __            |  |
//                |      |  |  PILLARS  |__|
//                |      |  |
//                |      |  |
//      INLET     |      |  |
//               \|/ h   |  |
//      _________________|  |___________________________


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    // The domain regions as defined in 'fsimicropillar.geo':
    int fluid = 1, solid = 2, inlet = 3, outlet = 4, sides = 5, clamp = 6;

    // The mesh can be curved!
    mesh mymesh("fsimicropillar.msh");

    // Define the fluid-solid interface:
    int fsinterface = regionintersection({fluid, solid});
    int boundary = regionunion({inlet,outlet,sides});

    // Confirm that the normal at the interface is pointing out of the pillar:
    normal(fsinterface).write(fsinterface, "normal.pos");


    // Field v is the flow velocity. It uses nodal shape functions "h1" with two components in 2D.
    // Field p is the relative pressure. Field u is the mechanical deflection. Field umesh stores
    // the smoothed mesh deformation while field y is the y coordinate.
    field v("h1xy"), p("h1"), u("h1xy"), umesh("h1xy"), y("y");

    // Force a no-slip (0 velocity) condition on the non-moving walls:
    v.setconstraint(sides);

    // Force a y-parabolic inflow velocity in the x direction increasing linearly over time at the inlet.
    // The channel height is h [m].
    double h = 120e-6;
    v.setconstraint(inlet, array2x1(   0.1 *4.0/(h*h)*y*(h-y) * t() , 0));
    // Set a 0 relative pressure [Pa] at the outlet:
    p.setconstraint(outlet);

    // The flag is clamped at the pillar bottom:
    u.setconstraint(clamp);


    // Use an order 1 interpolation for p and 2 for v on the fluid region (satisfies the BB condition).
    // Use order 2 for u on the solid region.
    p.setorder(fluid, 1); v.setorder(fluid, 2); u.setorder(solid, 2);
    umesh.setorder(fluid, 2); umesh.setorder(solid, 2);


    // Mesh deformation field umesh is forced to 0 on the fluid boundary, 
    // to u on the solid and smoothed with a Laplace formulation in the fluid.
    umesh.setconstraint(boundary); umesh.setconstraint(solid, u);

    // Classical Laplace formulation for each component:
    formulation laplacian;

    laplacian += integral(fluid, grad(dof(compx(umesh)))*grad(tf(compx(umesh))) + grad(dof(compy(umesh)))*grad(tf(compy(umesh))) );
    // This is needed to add the degrees of freedom in the solid region to the formulation:
    laplacian += integral(solid, dof(umesh)*tf(umesh) );


    // Dynamic viscosity of water [Pa.s] and density [kg/m3]:
    double mu = 8.9e-4, rhof = 1000;

    // Mechanical properties. Young's modulus E [Pa], Poisson's ratio nu [] and the density rho [kg/m3]:
    double E = 1e6, nu = 0.3, rhos = 1420;


    // Define the weak formulation for the fluid-structure interaction:
    formulation fsi;

    // No-slip condition. Force the fluid flow at the interface to dt(u):
    v.setconstraint(fsinterface, dt(u));
    // Add the force term applied by the fluid flow on the pillar (minus sign needed because the normal points outwards to the pillar).
    // Argument 'umesh' means the term is calculated on the mesh deformed by umesh.
    fsi += integral(fsinterface, umesh, -normal(fsinterface) * p * tf(u) );
    
    // The fluid velocity gradient has to be calculated on the 2D fluid elements first before being added as load to the interface elements.
    expression viscousforcetensor = mu*( grad(v)+transpose(grad(v)) );
    // Project each of its rows on a field to store the evaluated value. Only define degrees of freedom at the interface.
    field row1("h1xy"), row2("h1xy");
    row1.setorder(fluid, 2); row2.setorder(fluid, 2);
    fsi += integral(fluid, umesh, dof(row1, fsinterface)*tf(row1, fsinterface) - compx(viscousforcetensor)*tf(row1, fsinterface));
    fsi += integral(fluid, umesh, dof(row2, fsinterface)*tf(row2, fsinterface) - compy(viscousforcetensor)*tf(row2, fsinterface));
    
    // Use the fields as the viscous force tensor calculated on the fluid elements:
    fsi += integral(fsinterface, umesh, array2x1(dof(row1)*normal(fsinterface), dof(row2)*normal(fsinterface)) * tf(u) );


    // Classical elasticity with small-strain geometric-nonlinearity (obtained with the extra u argument). Update argument 0.0 for prestress.
    fsi += integral(solid, predefinedelasticity(dof(u), tf(u), u, E, nu, 0.0, "planestrain"));
    // Add the mechanic inertia term:
    fsi += integral(solid, -rhos*dtdt(dof(u))*tf(u));

    // Define the weak formulation for time-dependent incompressible laminar flow:
    fsi += integral(fluid, umesh, predefinednavierstokes(dof(v), tf(v), v, dof(p), tf(p), mu, rhof, 0, 0, true) );


    // Define the object for a generalized alpha time resolution (default parameters leads to Newmark).
    // An all zero initial guess for the fields and their time derivatives is set with 'vec(fsi)'.
    genalpha ga(fsi, vec(fsi), vec(fsi), vec(fsi));
    // Set the relative tolerance on the inner nonlinear iteration:
    ga.settolerance(1e-4);

    // A Laplace resolution is solved in the inner nonlinear loop of the gen alpha time resolution for mesh smoothing:
    ga.postsolve({laplacian});

    // Run in time from 0 to 1 sec by steps of 200 ms:
    std::vector<vec> sols = ga.runnonlinear(0, 0.1, 0.3)[0];

    for (int i = 0; i < sols.size(); i++)
    {
        // Transfer the solution at the ith timestep to the fields:
        v.setdata(fluid, sols[i]);
        u.setdata(solid, sols[i]);

        // Saving the velocity field on the deformed geometry requires the smoothed umesh field. 
        // It must be recalculated based on the current u (made available after the u.setdata call).
        solve(laplacian);

        // Write the fields to ParaView .vtk format:
        u.write(solid, "u" + std::to_string(1000 + i) + ".vtk", 2);
        // Write v on the deformed mesh:
        v.write(fluid, umesh, "v" + std::to_string(1000 + i) + ".vtk", 2);
    }

    // Output the max pillar deflection:
    double umax = norm(u).max(solid, 5)[0];
    std::cout << "Max pillar deflection: " << umax*1e6 << " um" << std::endl;

    // Code validation line. Can be removed.
    std::cout << (umax < 8.7115e-06 && umax > 8.7113e-06);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

