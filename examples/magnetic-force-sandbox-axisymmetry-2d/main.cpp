// This is a sandbox code to illustrate how magnetostatic stresses can be applied as a
// mechanical loading as well as how total magnetic forces can be computed. The geometry
// does not include corners to avoid the effect of singularities on the force calculation.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    setaxisymmetry();

    // Line region 'bnd' is defined below to perform a line integral of the Maxwell stress tensor.
    int bar = 1, ring = 2, air = 3, clamp = 4, sides = 5, bnd = 6;

    mesh mymesh;
    // Grow 5 element layers from the bar region:
    mymesh.selectlayer(7, air, bar, 5);
    // Take the layer skin:
    mymesh.selectskin(8, 7);
    // Take the bar skin:
    mymesh.selectskin(9, bar);
    // Remove the bar skin from the layer skin:
    mymesh.selectexclusion(bnd, 8, {9});
    // Refine by splitting once:
    mymesh.split(1);
    // Load the mesh:
    mymesh.load("ring2d.msh");

    int all = selectall();

    // Visualize the Maxwell stress tensor line integration region:
    expression(1).write(bnd, "bnd.vtu", 1);

    // Nodal shape functions 'h1' for the z component of the vector potential.
    field az("h1");

    // In 2D axisymmetry the vector potential only has a z component:
    expression a = array3x1(0,0,az);

    // Use interpolation order 4:
    az.setorder(all, 4);

    // Put a magnetic wall at the outer air sides:
    az.setconstraint(sides);
   
    // Vacuum magnetic permeability [H/m]: 
    double mu0 = 4.0*getpi()*1e-7;

    // Define the permeability in all regions.
    parameter mu;

    mu|all = mu0;
    // Overwrite on the magnetic bar region:
    mu|bar = 1e5*mu0;


    formulation magnetostatics;

    // The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
    magnetostatics += integral(all, 1/mu* curl(dof(a)) * curl(tf(a)) );
    // External current density source:
    magnetostatics += integral(ring, -array3x1(0,0,1e5) * tf(a));

    magnetostatics.solve();

    expression b = curl(a);

    b.write(all, "b.vtu", 2);
    az.write(all, "az.vtu", 2);
    
    // Compute the total force acting on the bar using a cell integral of the virtual work term:
    double fvolcalc = gettotalforce(bar, b/mu, mu)[1];
    // Compute the total force acting on the bar using a boundary integral of the Maxwell stress tensor T:
    expression T = 1/mu * ( b*transpose(b) - 0.5*b*b * eye(3) );
    double fsurcalc = 2*getpi()*compy(on(all, T)*normal(7)).integrate(bnd, 5);
    
    
    // Mechanical displacement field with 3 components:
    field u("h1xyz");
    
    // Use interpolation order 4:
    u.setorder(all, 4);
    
    // Clamp the x and z displacement components.
    // The y component is clamped using a port to be able to calculate the reaction force.
    u.compx().setconstraint(clamp);
    u.compz().setconstraint(clamp);
    
    // Define the y component displacement-force port pair:
    port Uy, Fy;
    u.compy().setport(clamp, Uy, Fy);

    formulation elasticity;
    
    // Force the y displacement on the clamp to zero:
    elasticity += Uy - 0.0;
    
    elasticity += integral(bar, predefinedelasticity(dof(u), tf(u), 100e9, 0.3));
    // Add the magnetostatic Maxwell stresses:
    elasticity += integral(all, predefinedmagnetostaticforce(tf(u, bar), b/mu, mu));
    // Alternatively an integration on the bar skin can be performed:
    // elasticity += integral(9, compx(on(air,T)) * normal(bar) * compx(tf(u)));
    // elasticity += integral(9, compx(on(bar,T)) * normal(air) * compx(tf(u)));
    // elasticity += integral(9, compy(on(air,T)) * normal(bar) * compy(tf(u)));
    // elasticity += integral(9, compy(on(bar,T)) * normal(air) * compy(tf(u)));
    
    elasticity.solve();
    
    u.write(bar, "u.vtu", 4);
    
    // Compute the total reaction force:
    double freaction = 2*getpi()*Fy.getvalue();

    std::cout << "Total force with cell integral is " << fvolcalc << " N" << std::endl;
    std::cout << "Total force with boundary integral is " << fsurcalc << " N" << std::endl;
    std::cout << "Mechanical reaction force at clamp is " << freaction << " N" << std::endl;
    
    double maxu = norm(u).max(bar, 5)[0];
    std::cout << "Max displacement is " << maxu << " m" << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (std::abs(maxu-4.3131e-8)/4.3131e-8 < 8e-5 && std::abs(fvolcalc-1496.3)/1496.3 < 2e-5 && std::abs(fsurcalc-1496.3)/1496.3 < 2e-5 && std::abs(freaction+1496.3)/1496.3 < 2e-5);
}

