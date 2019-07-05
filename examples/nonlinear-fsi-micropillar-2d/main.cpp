// This code simulates the fluid-structure interaction between an incompressible water 
// flow in a microchannel and two micropillars. The micropillars are modeled as elastic 
// structures. Small-strain geometric nonlinearity is taken into account.
// A monolithic fluid-structure coupling is used. A smooth mesh deformation is obtained
// by solving a Laplace formulation. 
// ALE can be used as an alternative to solve this FSI problem.
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

    // Load the mesh file (the mesh can be curved):
    mesh mymesh("fsimicropillar.msh");

    // Define the fluid-structure interface:
    int fsinterface = regionintersection({fluid, solid});

    // Confirm that the normal at the interface is pointing outwards of the pillar:
    normal(fsinterface).write(fsinterface, "normal.pos");


    // Field v is the flow velocity. It uses nodal shape functions "h1" with two components in 2D.
    // Field p is the relative pressure. Field u is the mechanical deflection. Field umesh stores
    // the smoothed mesh deformation while field y is the y coordinate.
    field v("h1xy"), p("h1"), u("h1xy"), umesh("h1xy"), y("y");

    // Force a no-slip (0 velocity) condition on the non-moving walls:
    v.setconstraint(sides);
    // Force the fluid flow at the fluid-structure interface to zero (no-slip condition for static simulation, dt(u) = 0):
    v.setconstraint(fsinterface);

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
    umesh.setconstraint(regionunion({inlet,outlet,sides})); umesh.setconstraint(solid, u);

    // Classical Laplace formulation for each component:
    formulation laplacian;

    laplacian += integral(fluid, grad(dof(compx(umesh)))*grad(tf(compx(umesh))) + grad(dof(compy(umesh)))*grad(tf(compy(umesh))) );
    // This is needed to add the degrees of freedom in the solid region to the formulation:
    laplacian += integral(solid, dof(umesh)*tf(umesh) );


    // Dynamic viscosity of water [Pa.s] and density [kg/m3]:
    double mu = 8.9e-4, rhof = 1000;

    // Mechanical properties. Young's modulus E [Pa], Poisson's ratio nu [] and the density rho [kg/m3]:
    double E = 1e6, nu = 0.3, rhos = 2000;

    
    // Calculate the viscous force tensor on the mesh deformed by field 'umesh'. The fluid velocity gradient 
    // has to be calculated on the 2D fluid elements first before being added as load to the fsi formulation.
    expression viscousforcetensor = mu*( grad(v) + transpose(grad(v)) );
    // Project each of the tensor entries to a field (use default order 1 interpolation, set to 2 for even more accuracy).
    field vft00("h1"), vft01("h1"), vft10("h1"), vft11("h1");
    
    formulation vftf00, vftf01, vftf10, vftf11;
    
    vftf00 += integral(fluid, umesh, dof(vft00)*tf(vft00) - entry(0,0,viscousforcetensor)*tf(vft00));
    vftf01 += integral(fluid, umesh, dof(vft01)*tf(vft01) - entry(0,1,viscousforcetensor)*tf(vft01));
    vftf10 += integral(fluid, umesh, dof(vft10)*tf(vft10) - entry(1,0,viscousforcetensor)*tf(vft10));
    vftf11 += integral(fluid, umesh, dof(vft11)*tf(vft11) - entry(1,1,viscousforcetensor)*tf(vft11));
    
    // Define the whole projected viscous force tensor for convenience:
    expression projectedviscousforcetensor = array2x2(vft00,vft01,vft10,vft11);
    

    // Define the weak formulation for the fluid-structure interaction:
    formulation fsi;
    
    // Classical elasticity with small-strain geometric nonlinearity (obtained with the extra u argument). Update argument 0.0 to add prestress.
    fsi += integral(solid, predefinedelasticity(dof(u), tf(u), u, E, nu, 0.0, "planestrain"));

    // Define the weak formulation for time-independent incompressible laminar flow (on the mesh deformed by 'umesh'):
    fsi += integral(fluid, umesh, predefinednavierstokes(dof(v), tf(v), v, dof(p), tf(p), mu, rhof, 0, 0, false) );

    // Add the hydrodynamic load pI - mu*( grad(v) + transpose(grad(v)) ) normal to the FSI interface (on the mesh deformed by 'umesh'):
    fsi += integral(fsinterface, umesh, ( p * -normal(fsinterface) - projectedviscousforcetensor * -normal(fsinterface) ) * tf(u) );


    double uprev, umax = 1;
    while (std::abs(umax-uprev)/std::abs(umax) > 1e-6)
    {
        // Project the viscous force tensor entries then solve the fluid-structure interaction formulation 
        // and the Laplace formulation to smooth the mesh. The fields are updated in each 'solve' call.
        solve({vftf00,vftf01,vftf10,vftf11});
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
    std::cout << (umax < 8.34643e-06 && umax > 8.34641e-06);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

