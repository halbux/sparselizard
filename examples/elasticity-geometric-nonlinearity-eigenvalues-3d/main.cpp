// This code gives the static deflection and mechanical vibration eigenfrequencies and eigenmodes 
// around that static deflection for a 3D bilayer disk (800um diameter, 10um thick) whose top layer is prestressed,
// whose bottom layer is clamped at its external boundary and on top of which the atmospheric pressure is pushing.
//
// Small-strain GEOMETRIC NONLINEARITY is used to simulate this prestressed membrane.
//
// With minimal extra effort the buckling of the membrane can be simulated with the genalpha time 
// resolution available in sparselizard.


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    // The domain regions as defined in 'disk.geo':
    int botlayer = 1, toplayer = 2, sur = 3, top = 4;
    
    // The mesh can be curved!
    mesh mymesh("disk.msh");
    
    // Define the whole domain for convenience:
    int vol = regionunion({botlayer, toplayer});
    
    // Nodal shape functions 'h1' with 3 components.
    // Field u is the membrane deflection.
    field u("h1xyz");
    
    // Use interpolation order 2 on the whole domain:
    u.setorder(vol, 2);
    
    // Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(sur);
    
    // E is Young's modulus [Pa]. nu is Poisson's ratio []. rho is the density [kg/m^3].
    // The material is a polymer.
    parameter E, nu, rho;
    E|botlayer = 10e9; nu|botlayer = 0.4; rho|botlayer = 1500;
    E|toplayer = 15e9; nu|toplayer = 0.4; rho|toplayer = 2000;  
    
    formulation elasticity;
    
    // The linear elasticity formulation with geometric nonlinearity for small strains is predefined:
    elasticity += integral(botlayer, predefinedelasticity(dof(u), tf(u), u, E, nu, 0.0), 0, 1);
    // The top layer is prestressed with 10 MPa in the x and y direction (sigma xx and yy):
    expression prestress(6,1,{10e6,10e6,0,0,0,0});
    elasticity += integral(toplayer, predefinedelasticity(dof(u), tf(u), u, E, nu, prestress), 0, 1);
    // Add the atmospheric pressure force at the top face (perpendicular to it).
    // With u provided as second argument all calculations are performed on the deformed 
    // mesh for this term (thus the normal is updated with the deflection).
    double pressure = 1e5;
    elasticity += integral(top, u, -pressure*normal(top)*tf(u), 0, 1);
    
    // Add the inertia terms to the formulation (required for eigenvalue computation).
    // The contribution below has a different tag (2) to be able to exclude it from 
    // the first .generate(1) call where the inertia terms should not be included.
    elasticity += integral(vol, -rho*dtdt(dof(u))*tf(u), 0, 2);
    
    
    ///// STEP 1: GET THE STATIC DEFLECTION WITH A NONLINEAR LOOP:
    
    // Start with an all-zero deflection:
    vec solu(elasticity);
    
    double relres = 1;
    while (relres > 1e-5)
    {
        // Generate and get the algebraic matrix A and vector b of the static Ax = b problem:
        elasticity.generate(1);
        
        mat A = elasticity.A();
        vec b = elasticity.b();
        
        // Compute the relative residual at the previous iteration:
        relres = (b-A*solu).norm() / b.norm();
        
        // Solve for x in Ax = b:
        solu = solve(A, b);
        // Make solution available in field u:
        u.setdata(vol, solu);
        
        // Write the deflection on the top surface of the membrane with an order 2 interpolation:
        u.write(top, "ustatic.pos", 2);
        
        std::cout << "Relative residual is: " << relres << std::endl;
    }
    
    std::cout << std::endl << "Peak static deflection: " << norm(u).max(vol,4)[0]*1e6 << " um" << std::endl;
    
    
    ///// STEP 2: COMPUTE THE EIGENFREQUENCIES AND EIGENMODES AROUND THE STATIC DEFLECTION:
    
    elasticity.generate({1,2});
    
    // Get the stiffness and mass matrix:
    mat K = elasticity.K();
    mat M = elasticity.M();
    
    // Remove the rows and columns corresponding to the 0 constraints:
    K.removeconstraints();
    M.removeconstraints();
    
    // Create the object to solve the generalized eigenvalue problem K*x = lambda*M*x :
    eigenvalue eig(K, M);
    
    // Compute the 5 eigenvalues closest to the target magnitude 0.0 (i.e. the 5 first ones):
    eig.compute(5, 0.0);
    
    // Print the eigenfrequencies:
    eig.printeigenfrequencies();
    
    // The eigenvectors are real thus we only need the real part:
    std::vector<vec> myeigenvectors = eig.geteigenvectorrealpart();
    
    // Loop on all eigenvectors found:
    for (int i = 0; i < myeigenvectors.size(); i++)
    {
        // Transfer the data from the ith eigenvector to field u:
        u.setdata(top, myeigenvectors[i]);
        // Write the deflection on the top surface of the membrane with an order 2 interpolation:
        u.write(top, "ueigenmode"+std::to_string(i)+".pos", 2);
    }

    // Code validation line. Can be removed.
    std::cout << (eig.geteigenvaluerealpart()[0] < 1.1453e+12 && eig.geteigenvaluerealpart()[0] > 1.1451e+12);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

