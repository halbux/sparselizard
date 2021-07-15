// This code simulates the static deflection of an electrically actuated 2D micromembrane.
// The problem is nonlinear because of the electrostatic-elastic coupling.
//
// The electrostatic problem and the electrostatic forces are computed on a mesh 
// deformed by the mechanical displacement. The mesh in the vacuum gap below the membrane 
// is deformed smoothly by solving an additional Laplace problem.
//
// The interpolation order for the electric potential field and mechanical displacement
// is adapted to every geometrical region.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'cmut2d.geo':
    int insulator = 1, pillars = 2, vacuumgap = 3, membrane = 4, ground = 5, clamp = 5, electrode = 6;
    
    // The mesh can be curved!
    mesh mymesh("cmut2d.msh");
    
    // Define the region for the mechanical problem:
    int solid = selectunion({insulator, pillars, membrane});
    // Define the region for the electric problem:
    int electricdomain = selectunion({insulator, pillars, membrane, vacuumgap});
    
    // Nodal shape functions 'h1' for the electric potential field v
    // and membrane deflection u (u has components x and y).
    // umesh smoothly deforms the mesh in the vacuum gap for a given membrane deflection.
    field u("h1xy"), umesh("h1xy"), v("h1");

    // Use interpolation order:
    //
    // - 3 for u on the membrane
    // - 2 for u on the pillars
    // - 1 elsewhere
    //
    // - 1 for the electric potential v
    
    u.setorder(membrane, 3);
    u.setorder(vacuumgap, 3);
    u.setorder(pillars, 2);
    u.setorder(insulator, 1);
    
    v.setorder(electricdomain, 1);
    
    umesh.setorder(solid, 1);
    umesh.setorder(electricdomain, 1);
    
    // Clamp and ground (i.e. 0 valued-Dirichlet conditions for u and v):
    u.setconstraint(clamp);
    v.setconstraint(ground);
    // Force the electric potential on the electrode to a close-to-pull-in voltage:
    v.setconstraint(electrode, 200);
  
    // E is Young's modulus. nu is Poisson's ratio.
    // epsilon is the electric permittivity. 
    //
    // The membrane is polysilicon, the insulator is silicon dioxyde.
    
    parameter E, nu, epsilon;
    
    E|insulator =		70e9;
    E|pillars =			150e9;
    E|membrane =		150e9;
    
    nu|insulator =		0.17;
    nu|pillars =		0.3;
    nu|membrane =		0.3;
    
    epsilon|vacuumgap =	8.854e-12;
    epsilon|insulator =	3.9*8.854e-12;
    epsilon|pillars =	11.7*8.854e-12;
    epsilon|membrane =	11.7*8.854e-12;
    
    // An electrostatic formulation is used for the electric problem.
    // An elasticity formulation is used for the mechanical problem.
    formulation electrostatics, elasticity;
	
    // Weak electrostatic formulation, computed on the mesh deformed by field umesh:
    electrostatics += integral(electricdomain, umesh, epsilon*grad(dof(v))*grad(tf(v)));
    
    // The linear elasticity formulation is classical and thus predefined:
    elasticity += integral(solid, predefinedelasticity(dof(u), tf(u), E, nu, "planestrain"));
    
    // Electrostatic forces, computed on the elements of the whole electric domain
    // but with mechanical deflection test functions tf(u) only on solid.
    // 
    // The electrostatic forces often appear in MEMS simulations and are thus predefined.
    // The inputs are the gradient of the test function of u defined on the mechanical domain,
    // the gradient of the previously computed electric potential field and the electric permittivity.
    //
    // The electrostatic forces are computed on the mesh deformed by field umesh.
    elasticity += integral(electricdomain, umesh, predefinedelectrostaticforce(tf(u,solid), grad(v), epsilon));
    
    
    // Solve the Laplace equation in the vacuum gap to smoothly deform the mesh.
    // umesh is forced to field u on region solid:
    umesh.setconstraint(solid, u);

    formulation laplacian;
    laplacian += integral(vacuumgap, grad(dof(compx(umesh)))*grad(tf(compx(umesh))) + grad(dof(compy(umesh)))*grad(tf(compy(umesh))) );
    
    
    // NONLINEAR ITERATION TO GET THE STATIC DEFLECTION:
  	
    // Start with an all-zero solution vector for the elasticity formulation:
    vec solu(elasticity);
    
    double relresnorm = 1; int iter = 0; 
    while (relresnorm > 1e-5)
    {
        electrostatics.generate();
        vec solv = solve(electrostatics.A(), electrostatics.b());
        // Transfer the data from the solution vector to the v field:
        v.setdata(electricdomain, solv);
        // Write the electric field with an order 1 interpolation.
        // The electric field is computed and saved on the geometry deformed by umesh.
        (-grad(v)).write(electricdomain, umesh, "E.pos", 1);
        
        // Use the now known electric potential v to compute the membrane deflection:
        elasticity.generate();
        
        vec b = elasticity.b();
        mat A = elasticity.A();
        // Compute the norm of the relative residual:
        relresnorm = (b-A*solu).norm()/b.norm();

        solu = solve(A,b);
        u.setdata(solid, solu);
        // Write the deflection u with an order 3 interpolation:
        u.write(solid, "u.pos", 3); 
        
        // Smooth the mesh on the vacuum gap:
        laplacian.generate();
        vec solumesh = solve(laplacian.A(), laplacian.b());
        // Save the smoothed mesh on the vacuum region:
        umesh.setdata(vacuumgap, solumesh);
        // Also save the u field on region solid to umesh.
        // This is done by selecting field u with |u on the solu vector.
        umesh.setdata(solid, solu|u);
        umesh.write(electricdomain, "umesh.pos", 3);
        
        // Print the iteration number and relative residual norm:
        std::cout << "Relative residual norm @iteration " << iter << " is " << relresnorm << std::endl;
        iter++;
    }
    
    // Code validation line. Can be removed.
    std::cout << (compy(grad(v)).integrate(vacuumgap, u, 4) < 0.0022905 && compy(grad(v)).integrate(vacuumgap, u, 4) > 0.0022901);
}

