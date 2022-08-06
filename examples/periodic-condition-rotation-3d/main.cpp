// This simple code shows how to setup a 3D mechanical problem with a periodic condition.
// A central 1.3 um thick, 50 um diameter polysilicon micromembrane stands above a 300 nm deep cavity.
// The central micromembrane is surrounded by 6 identical ones, one every 60 degrees.
// Assuming a periodic vibration with 60 degrees periodicity allows to reduce the problem to only a 1/6 of the geometry.
// 
// In the .nas mesh used the meshes on both faces of the periodic condition do not match.


#include "sparselizard.h"


using namespace sl;

void processmesh(void);

int main(void)
{	
    // Name the regions for the inner and outer electrode, the clamp and the regions 'gamma' on which to apply the periodic condition:
    int electrodein = 1, electrodeout = 2, clamp = 3, gamma1 = 4, gamma2 = 5, cavity = 4007, solid = 4008, all = 4009;

    // Define all regions needed in the source .nas mesh and save it in .msh format.
    processmesh();

    mesh mymesh("cmutperiodic.msh");

    // The periodic condition is only applied to the solid region:
    gamma1 = selectintersection({gamma1, solid}, 2);
    gamma2 = selectintersection({gamma2, solid}, 2);


    wallclock clk;

    // Harmonic simulation at f0 = 1 MHz:
    setfundamentalfrequency(1e6);

    // Nodal shape functions 'h1' with 3 components for u, the membrane deflection.
    // Use harmonic 2 to have u = U*sin(2pi*f0*t).
    field u("h1xyz", {2});

    // Use interpolation order 2 everywhere:
    u.setorder(all, 2);

    // Clamp on surface 'clamp' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);

    // E is Young's modulus. nu is Poisson's ratio. 
    double E = 160e9, nu = 0.22, rho = 2320;


    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(solid, predefinedelasticity(dof(u), tf(u), E, nu));
    // Add a pressure load at frequency f0 on both inner and outer electrodes:
    double p = 1e5;
    elasticity += integral(electrodein, -p*compz(tf(u.harmonic(2))));
    elasticity += integral(electrodeout, -p*compz(tf(u.harmonic(2))));
    // Add the inertia term:
    elasticity += integral(solid, -rho*dtdt(dof(u))*tf(u));


    // Add the periodic condition between gamma1 and gamma2.
    // Region gamma2 is obtained from gamma1 by a 60 degrees rotation around z (rotation center is the origin).
    elasticity += periodicitycondition(gamma1, gamma2, u, {0,0,0}, {0,0,60}, 1.0, 1);


    // Generate, solve and store the solution to field u:
    elasticity.solve();

    // Write the deflection to ParaView .vtk format. Write with an order 2 interpolation.
    u.write(solid, "u.vtk", 2);

    // Confirm that the periodic condition is correct by comparing the inner and outer cavity deflection:
    double ucenterin = 1e9 * norm(u.harmonic(2)).interpolate(solid,{0,0,1.5e-6})[0];
    double ucenterout = 1e9 * norm(u.harmonic(2)).interpolate(solid,{60e-6,0,1.5e-6})[0];

    std::cout << "Deflection at center of inner/outer cavity is " << ucenterin << " / " << ucenterout << " nm" << std::endl;

    clk.print();

    // Code validation line. Can be removed.
    std::cout << (ucenterin < 26.5977 && ucenterin > 26.5975);
}

void processmesh(void)
{
    // Define the central electrode, outer electrode and clamp regions as well as the regions to apply the periodic condition.
    int elecc = 1, eleco = 2, clamp = 3, gamma1 = 4, gamma2 = 5;

    setphysicalregionshift(1000);

    mesh mymesh1;

    mymesh1.load("cmutperiodic.nas", 0);

    int vac = selectunion({4001,4005});
    int solid = selectunion({4002,4003,4004,4006});
    int all = selectunion({vac,solid});

    // Rotate the mesh to easily select the bottom side for the periodic condition:
    mymesh1.rotate(0,0,30);
    mymesh1.write("cmutperiodic.msh");
    
    setphysicalregionshift(0);

    mesh mymesh2;

    mymesh2.selectbox(elecc, 4001, 2, {-10,10,-10,10,0.3e-6,0.3e-6});
    mymesh2.selectbox(eleco, 4006, 2, {-10,10,-10,10,0.3e-6,0.3e-6});
    mymesh2.selectbox(clamp, all, 2, {-10,10,-10,10,0,0});
    mymesh2.selectbox(gamma1, all, 2, {-10,10,0,0,-10,10});

    mymesh2.load("cmutperiodic.msh", 0);

    // Rotate to the other direction to align the other region for the periodic condition:
    mymesh2.rotate(0,0,-60);
    mymesh2.write("cmutperiodic.msh");

    mesh mymesh3;

    mymesh3.selectbox(gamma2, all, 2, {-10,10,0,0,-10,10});

    mymesh3.load("cmutperiodic.msh", 0);

    // Bring the mesh back to its original angle:
    mymesh3.rotate(0,0,30);

    // Write the processed mesh:
    mymesh3.write("cmutperiodic.msh");
}

