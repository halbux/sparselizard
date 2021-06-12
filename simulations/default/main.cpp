#include "sparselizard.h"

using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'disk.geo':
    int vol = 1, sur = 2, top = 3;
    
    // The mesh can be curved!
    mesh mymesh("disk.msh");
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the membrane deflection.
    field u("h1xyz");

    // Use interpolation order 2 on 'vol', the whole domain:
    u.setorder(vol, 2);
    
    // Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(sur);
  
    // E is Young's modulus. nu is Poisson's ratio. 
    parameter E, nu;
    E|vol = 150e9; nu|vol = 0.3;
  
    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), E, nu));
    // Add a volumic force in the -z direction:
    elasticity += integral(vol, array3x1(0,0,-10)*tf(u));

    // Generate, solve and transfer the solution to field u:
    elasticity.solve();
    
    // Write the deflection to ParaView .vtk format.
    // Write with an order 2 interpolation. Exaggerate the deflection by a factor 0.5e9.
    (0.5e9*u).write(vol, "u.vtk", 2);
}

