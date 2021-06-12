// This example shows how to load a mesh from the GMSH API.
// The GMSH API must be available (run the appropriate script in 'install_external_libs' for that).


#include "sparselizard.h"
#include "gmsh.h"


using namespace sl;

int main(void)
{	
    // EXAMPLE 1.
    // Load from file via GMSH API and solve the elasticity equations.
    
    gmsh::initialize();
    
    // The domain regions as defined in 'disk.geo':
    int vol = 1, clamp = 2;
    
    gmsh::open("disk.msh");
    // Load mesh available in GMSH API.
    // Set verbosity to 2 to print the physical regions info.
    mesh mymesh("gmsh:api", 2);
    
    gmsh::finalize();


    // EXAMPLE 2. BASED ON A GEOMETRY PROPOSED IN THE GMSH EXAMPLES.
    // Create mesh with GMSH API, load it and solve the elasticity equations.
    /*
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);

    gmsh::model::add("boolean");

    gmsh::option::setNumber("Mesh.Algorithm", 6);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0.4);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.4);

    double R = 1.4, Rs = R*.7, Rt = R*1.25;

    std::vector<std::pair<int, int> > ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;
    gmsh::model::occ::addBox(-R,-R,-R, 2*R,2*R,2*R, 1);
    gmsh::model::occ::addSphere(0,0,0,Rt, 2);
    gmsh::model::occ::intersect({{3, 1}}, {{3, 2}}, ov, ovv, 3);
    gmsh::model::occ::addCylinder(-2*R,0,0, 4*R,0,0, Rs, 4);
    gmsh::model::occ::addCylinder(0,-2*R,0, 0,4*R,0, Rs, 5);
    gmsh::model::occ::addCylinder(0,0,-2*R, 0,0,4*R, Rs, 6);
    gmsh::model::occ::fuse({{3, 4}, {3, 5}}, {{3, 6}}, ov, ovv, 7);
    gmsh::model::occ::cut({{3, 3}}, {{3, 7}}, ov, ovv, 8);

    gmsh::model::occ::synchronize();
    
    // Add new physical region numbers.
    int vol = gmsh::model::addPhysicalGroup(3, {8}, -1);
    int clamp = gmsh::model::addPhysicalGroup(2, {6}, -1);

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::refine();
    gmsh::model::mesh::setOrder(2);
    //gmsh::model::mesh::partition(4);

    gmsh::write("boolean.msh");

    // Load mesh from GMSH API to sparselizard.
    mesh mymesh("gmsh:api");

    gmsh::finalize();
    */
    
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the mechanical deflection.
    field u("h1xyz");

    // Use interpolation order 2 on 'vol', the whole domain:
    u.setorder(vol, 2);
    
    // Clamp on surface 'clamp' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);
  
    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), 150e9, 0.3));
    // Add a volumic force in the -z direction:
    elasticity += integral(vol, array1x3(0,0,-10)*tf(u));

    elasticity.generate();

    vec solu = solve(elasticity.A(), elasticity.b());

    // Transfer the data from the solution vector to the u field:
    u.setdata(vol, solu);
    // Write the deflection to ParaView .vtk format. Write with an order 2 interpolation.
    u.write(vol, "u.vtk", 2);
    
    // Code validation line. Can be removed.
    std::cout << (compz(u).integrate(vol, u, 5) < -1.24776e-10 && compz(u).integrate(vol, u, 5) > -1.24780e-10);
}

