// This code shows the most accurate way to compute mechanical reaction forces.
// Here the reaction forces of a thin 3D disk are calculated at the clamping.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'disk.geo':
    int vol = 1, sur = 2, top = 3;
    
    // The mesh is curved to capture accurately the circular geometry:
    mesh mymesh("disk.msh");
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the membrane deflection.
    field u("h1xyz");
    
    // Volume of the disk (radius is 1 m):
    double volume = getpi()*0.1;

    // Use interpolation order 2 on 'vol', the whole domain:
    u.setorder(vol, 2);
    
    // Associate a U/F (displacement/force) port pair to field u on the clamp face 'sur'.
    // The displacement components on the port region must be constant (here equal to zero).
    port Ux, Uy, Uz, Fx, Fy, Fz;
    u.compx().setport(sur, Ux, Fx);
    u.compy().setport(sur, Uy, Fy);
    u.compz().setport(sur, Uz, Fz);
    //                      |  |
    //            primal port  dual port
    //
    // The dual port holds the global Neumann term on the port region.
    // For an elasticity formulation this equals the reaction force.
  
    formulation elasticity;

    // Add port relations Ux = 0, Uy = 0 and Uz = 0 to clamp the disk:
    elasticity += Ux - 0.0;
    elasticity += Uy - 0.0;
    elasticity += Uz - 0.0;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), 150e9, 0.3));
    // Add an arbitrary uniform volumic force (N/m^3) in the x, y and z direction:
    elasticity += integral(vol, array3x1(1,2,3)*tf(u));

    // Generate, solve and transfer the solution to field u and to the ports:
    elasticity.solve();

    std::cout << "Reaction force on clamp is (" << Fx.getvalue() << ", " << Fy.getvalue() << ", " << Fz.getvalue() << ") N" << std::endl;
    std::cout << "Exact reaction on clamp is (" << -1*volume << ", " << -2*volume << ", " << -3*volume << ") N" << std::endl;
    
    // Print the max deflection:
    std::cout << "Max displacement is " << norm(u).max(vol, 5)[0] << " m" << std::endl;
    
    // Write the deflection to ParaView .vtk format. Write with an order 2 interpolation.
    u.write(vol, "u.vtk", 2);
    
    // Code validation line. Can be removed.
    double fsum = Fx.getvalue() + Fy.getvalue() + Fz.getvalue();
    std::cout << (-fsum < 1.8849559 && -fsum > 1.8849555);
}

