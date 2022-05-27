// This code gives the mechanical deflection of a 2D truss with a pressure load applied on one of its sides. 
// The pressure direction is updated as the truss deforms. Plane stress is assumed.
// Geometric nonlinearity is taken into account.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'truss2d.geo':
    int solid = 1, clamp = 2, load = 3;

    // Load the GMSH 4 format mesh with gmsh:
    mesh mymesh("gmsh:truss2d.msh", 1);

    // Nodal shape functions 'h1' with 2 components for the mechanical displacement:
    field u("h1xy");

    // Use interpolation order 2 on the whole domain:
    u.setorder(solid, 2);

    // Clamp the truss (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);

    // E [Pa] is Young's modulus. nu [] is Poisson's ratio. 
    double E = 1e9, nu = 0.3;

    formulation elasticity;

    // The elasticity formulation with geometric nonlinearity is classical and thus predefined:
    elasticity += integral(solid, predefinedelasticity(dof(u), tf(u), u, E, nu, 0.0, "planestress"));
    // Add a pressure force on the 'load' line. Compute the force on the mesh deformed by field u.
    // The normal direction moves with the mesh due to the 'u' argument.
    elasticity += integral(load, u, -0.8e3*normal(solid)*tf(u));

    double prevumax = 1, umax = 2;
    while (std::abs(umax-prevumax)/std::abs(prevumax) > 1e-8)
    {
        elasticity.solve();

        prevumax = umax;
        umax = norm(u).max(solid, 5)[0];

        std::cout << "Max u: " << umax << " m (rel change " << std::abs(umax-prevumax)/std::abs(prevumax) << ")" << std::endl;
    }

    // Write the deflection to ParaView .vtk format with an order 2 interpolation.
    u.write(solid, "u.vtk", 2);

    // Code validation line. Can be removed.
    std::cout << (umax < 0.0409512 && umax > 0.0409510);
}

