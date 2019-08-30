// This code gives the mechanical vibration eigenfrequencies and eigenmodes of a 
// 3D disk that is clamped at its outer face and subject to proportional damping.
// A second order polynomial eigenvalue problem is solved.


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
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

    // E is Young's modulus. nu is Poisson's ratio. rho is the volumic mass.
    parameter E, nu, rho;
    E|vol = 150e9; nu|vol = 0.3; rho|vol = 2330;

    formulation elasticity;

    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(vol, predefinedelasticity(dof(u), tf(u), E, nu));
    // Add the inertia terms:
    elasticity += integral(vol, -rho*dtdt(dof(u))*tf(u));

    elasticity.generate();

    // Get the stiffness and mass matrix:
    mat K = elasticity.K();
    mat M = elasticity.M();

    // Remove the rows and columns corresponding to the 0 constraints:
    K.removeconstraints();
    M.removeconstraints();

    // Proportional damping:
    double alpha = 0.00003, beta = 0.00006;
    mat C = alpha*M + beta*K;

    // Create the object to solve the second order polynomial eigenvalue problem (M*lambda^2 + C*lamda + K)x = 0 :
    eigenvalue eig(K, C, M); // Replace by eig(K, M) for an undamped simulation

    // Compute the 10 eigenvalues closest to the target magnitude 0.0 (i.e. the 10 first ones):
    eig.compute(10, 0.0);

    // Print the eigenfrequencies:
    eig.printeigenfrequencies();

    // Create a harmonic field to store the real and imaginary eigenvector part:
    field ueig("h1xyz", {2,3});
    ueig.setorder(vol, 2);

    // The eigenvectors are real thus we only need the real part:
    std::vector<vec> myrealeigenvectors = eig.geteigenvectorrealpart();
    std::vector<vec> myimageigenvectors = eig.geteigenvectorimaginarypart();

    // Loop on all eigenvectors found:
    for (int i = 0; i < myrealeigenvectors.size(); i++)
    {
        // Transfer the data from the ith eigenvector to field ueig.
        // Use |u to select field u in the vector structure since ueig is not part of it:
        ueig.harmonic(2).setdata(vol, myrealeigenvectors[i]|u);
        ueig.harmonic(3).setdata(vol, myimageigenvectors[i]|u);
        // Write the deflection on the top surface of the membrane with an order 2 interpolation:
        ueig.write(vol, "ueig"+std::to_string(i)+".vtk", 2);
    }

    // Code validation line. Can be removed.
    double checkval = eig.geteigenvaluerealpart()[1]*myimageigenvectors[3].norm();
    std::cout << (checkval > -0.0021397 && checkval < -0.0021396);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

