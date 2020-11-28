// This code shows how to load a mesh file to shape objects for editing.
// The 'disk.msh' file is loaded and a thin extra disk slice is added on top of it.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'disk.msh':
    int vol = 1, sur = 2, top = 3;
    
    // Load the mesh 'disk.msh' for editing.
    // 'diskshapes[d]' lists the shapes containing every physical 
    // region of dimension d (0D, 1D, 2D or 3D) in file 'disk.msh'.
    std::vector<std::vector<shape>> diskshapes = loadshape("disk.msh");
    
    // Add a thin slice on top of the disk (diskshapes[2][1] is the disk top face):
    shape thinslice = diskshapes[2][1].extrude(vol, 0.02, 2);
    
    // Load all shapes of interest (the disk and the additional thin slice):
    mesh mymesh({diskshapes[2][0], diskshapes[2][1], diskshapes[3][0], thinslice});
    mymesh.write("editeddisk.msh");
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the membrane deflection.
    field u("h1xyz");

    // Use interpolation order 2 on 'vol', the whole domain:
    u.setorder(vol, 2);
    
    // Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(sur);
  
    // E is Young's modulus. nu is Poisson's ratio. 
    double E = 150e9, nu = 0.3;
  
    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), E, nu));
    // Add a volumic force in the -z direction:
    elasticity += integral(vol, array1x3(0,0,-10)*tf(u));

    elasticity.generate();

    vec solu = solve(elasticity.A(), elasticity.b());

    // Transfer the data from the solution vector to the u field:
    u.setdata(vol, solu);
    // Write the deflection to ParaView .vtk format.
    // Write with an order 2 interpolation. Exaggerate the deflection by a factor 0.5e9.
    (0.5e9*u).write(vol, "u.vtk", 2);
    
    // Code validation line. Can be removed.
    double maxu = norm(u).max(vol, 5)[0];
    std::cout << (maxu < 8.69951e-10 && maxu > 8.69947e-10);
}

