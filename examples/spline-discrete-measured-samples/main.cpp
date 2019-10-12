#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    // The domain regions as defined in 'disk.geo':
    int vol = 1, sur = 2, top = 3;
    
    // The mesh can be curved!
    mesh mymesh("disk.msh");
    
    // Nodal shape functions 'h1' with 3 components for the mechanical displacement u [m]. 
    // Field T is the temperature [K] and x/y is the x/y coordinate.
    field u("h1xyz"), T("h1"), x("x"), y("y");

    // Use interpolation order 2 on 'vol', the whole domain:
    u.setorder(vol, 2);
    
    // Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(sur);
  
    // Load the measured Young's modulus versus temperature data samples in a spline object:
    spline measureddata("steel-stiffness-temperature.txt");
    // Define the expression giving Young's modulus [Pa] as a spline interpolation of the temperature data:
    expression E(measureddata, T);

    // nu is Poisson's ratio []. 
    double nu = 0.3;
    
    // Define am arbitrary space-dependent temperature field for illustration:
    T.setvalue(vol, 473+100*(1+x)*(1+y));
    
  
    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), E, nu));
    // Add a volumic force in the -z direction:
    elasticity += integral(vol, array1x3(0,0,-10)*tf(u));

    elasticity.generate();

    vec solu = solve(elasticity.A(), elasticity.b());

    // Transfer the data from the solution vector to the u field:
    u.setdata(vol, solu);
    // Write the deflection to ParaView .vtk format with an order 2 interpolation:
    u.write(vol, "u.vtk", 2);
    // Write Young's modulus in space for illustration:
    E.write(vol, "E.vtk", 2);

    // Print the peak deflection:
    double umax = norm(u).max(vol,5)[0];
    std::cout << umax << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (umax < 9.63876e-10 && umax > 9.63874e-10);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

