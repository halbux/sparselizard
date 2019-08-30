// This code gives the mechanical vibration eigenfrequencies and eigenmodes of a 
// 3D disk that is clamped at its outer face and subject to proportional damping.


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
    // This defines the damping matrix C. It can alternatively be obtained with elasticity.C()
    // after having defined an additional damping term (i.e. a term with a dt(dof(u))) in the formulation.
    mat C = alpha*M + beta*K;

    // Create the object to solve the second order polynomial eigenvalue problem (M*lambda^2 + C*lambda + K)u = 0 :
    eigenvalue eig(K, C, M); // Replace by eig(K, M) for an undamped simulation

    // Compute the 10 eigenvalues closest to the target magnitude 0.0 (i.e. the 10 first ones):
    eig.compute(10, 0.0);

    // Print the eigenfrequencies:
    eig.printeigenfrequencies();

    // The eigenvectors are real only in the undamped case:
    std::vector<vec> myrealeigenvectors = eig.geteigenvectorrealpart();
    std::vector<vec> myimageigenvectors = eig.geteigenvectorimaginarypart();

    // Loop on all eigenvectors found:
    for (int i = 0; i < myrealeigenvectors.size(); i++)
    {
        // Transfer the data from the ith eigenvector to field u:
        u.setdata(vol, myrealeigenvectors[i]);
        // Write the deflection on the membrane with an order 2 interpolation:
        u.write(vol, "u"+std::to_string(i)+".vtk", 2);
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

