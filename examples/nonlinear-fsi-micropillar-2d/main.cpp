// WARNING: THIS RESOLUTION INCLUDES 3X MORE UNKNOWNS THAN NEEDED AND RUNS SLOW FOR NOW
// A FUNCTION TO AVOID THE NEED OF THESE 3X UNKNOWNS WILL SOON BE ADDED!!!!!!!!!!

// This code simulates the fluid-structure interaction between an incompressible water 
// flow in a microchannel and two micropillars. The micropillars are modeled as elastic 
// structures. Small-strain geometric nonlinearity is taken into account.
// A monolithic fluid-structure coupling is used. A smooth mesh deformation is obtained
// by solving a Laplace formulation. Alternatively an ALE resolution can be used.
//
// A parabolic normal flow velocity is forced at the inlet and 
// a zero pressure is imposed at the outlet.
//
//
// The microchannel geometry is as follows (l = 350 um, h = 120 um):
//
//                                       l
//      <------------------------------------------------------------------>
//      ________________________________    ________________________________
//                                      |  |
//               /|\                    |  |                          OUTLET
//                |                     |  |
//                |                     |  |
//                |       __            |  |
//                |      |  |  PILLARS  |__|
//                |      |  |
//                |      |  |
//      INLET     |      |  |
//               \|/ h   |  |
//      _________________|  |_______________________________________________


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

    // Force a y-parabolic inflow velocity in the x direction at the inlet.
    // The channel height is h [m].
    double h = 120e-6;
    v.setconstraint(inlet, array2x1( 0.03 * 4.0/(h*h)*y*(h-y), 0));
    // Set a 0 relative pressure [Pa] at the outlet:
    p.setconstraint(outlet);

    // The pillars are clamped at their bottom side:
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
    double E = 1e6, nu = 0.3, rhos = 2000;


    // Define the weak formulation for the fluid-structure interaction:
    formulation fsi;

    // No-slip condition. Force the fluid flow at the fluid-structure interface to zero:
    v.setconstraint(fsinterface);
    // Add the force term applied by the fluid flow on the pillar (minus sign needed because the normal points outwards to the pillars).
    // Argument 'umesh' means the term is calculated on the mesh deformed by umesh.
    fsi += integral(fsinterface, umesh, -normal(fsinterface) * p * tf(u) );
    
    // The fluid velocity gradient has to be calculated on the 2D fluid elements first before being added as load to the interface elements.
    expression viscousforcetensor = mu*( grad(v)+transpose(grad(v)) );
    // Project each of its rows on a field to store the evaluated value.
    field row1("h1xy"), row2("h1xy");
    row1.setorder(fluid, 2); row2.setorder(fluid, 2);
    fsi += integral(fluid, umesh, dof(row1)*tf(row1) - compx(viscousforcetensor)*tf(row1));
    fsi += integral(fluid, umesh, dof(row2)*tf(row2) - compy(viscousforcetensor)*tf(row2));
    
    // Use the fields as the viscous force tensor calculated on the fluid elements:
    fsi += integral(fsinterface, umesh, array2x1(dof(row1)*normal(fsinterface), dof(row2)*normal(fsinterface)) * tf(u) );


    // Classical elasticity with small-strain geometric nonlinearity (obtained with the extra u argument). Update argument 0.0 for prestress.
    fsi += integral(solid, predefinedelasticity(dof(u), tf(u), u, E, nu, 0.0, "planestrain"));

    // Define the weak formulation for time-independent incompressible laminar flow:
    fsi += integral(fluid, umesh, predefinednavierstokes(dof(v), tf(v), v, dof(p), tf(p), mu, rhof, 0, 0, false) );


    double uprev, umax = 1;
    while (std::abs(umax-uprev)/std::abs(umax) > 1e-6)
    {
        // Solve the fluid-structure interaction formulation and the Laplace formulation to smooth the mesh. 
        // The fields are updated in the 'solve' call.
        solve(fsi);
        solve(laplacian);
        
        // Output the max pillar deflection:
        uprev = umax; 
        umax = norm(u).max(solid, 5)[0];
        std::cout << "Max pillar deflection (relative change): " << umax*1e6 << " um (" << std::abs(umax-uprev)/std::abs(umax) << ")" << std::endl;
    }
    
    // Write the fields to ParaView .vtk format:
    u.write(solid, "u.vtk", 2);
    // Write v on the deformed mesh:
    v.write(fluid, umesh, "v.vtk", 2);

    // Code validation line. Can be removed.
    std::cout << (umax < 8.35682e-06 && umax > 8.35680e-06);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

