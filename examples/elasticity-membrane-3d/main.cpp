// This code simulates the static mechanical deflection of a 3D disk when a volume
// force is applied. The disk is clamped at its outer face.


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    // The domain regions as defined in 'circle.geo':
    int vol = 1, sur = 2, top = 3;
    
    // The mesh can be curved!
    mesh mymesh("disk.msh");
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the membrane deflection.
    field u("h1xyz");

    // Use interpolation order 3 on 'vol', the whole domain:
    u.setorder(vol, 3);
    
    // Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(sur);
  
    // E is Young's modulus. nu is Poisson's ratio. 
    parameter E, nu;
    E|vol = 150e9; nu|vol = 0.3;
  
    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(u, E, nu));
    // Add a volumic force in the -z direction:
    elasticity += integral(vol, array1x3(0,0,-10)*tf(u));

    elasticity.generate();

    vec solu = solve(elasticity.A(), elasticity.b());

    // Transfer the data from the solution vector to the u field:
    u.getdata(vol, solu);
    // Write the deflection on the top surface of the membrane.
    // Write with an order 3 interpolation. Exagerate the deflection by a factor 1e9.
    (0.5e9*u).write(top, "u.pos", 3);
    
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}









