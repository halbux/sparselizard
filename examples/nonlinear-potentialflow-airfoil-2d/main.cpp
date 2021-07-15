// This code simulates a potential flow (subsonic flow) around a horizontal NACA0012 airfoil.
// The problem is nonlinear because the air density depends on the air speed.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'airfoil2D.geo':
    int air = 1, airfoil = 2, downstream = 3, upstream = 4;
    
    // Load the 'airfoil2D' mesh: 
    mesh mymesh("airfoil2D.msh");
    
    // Define the whole outer boundary:
    int gammaouter = selectunion({upstream, downstream});
    
    // Define the velocity potential field 'phi' with standard nodal shape functions ("h1").
    // grad(phi) is the fluid velocity. Field x is the x coordinate field.
    field phi("h1"), x("x");
    
    // Use interpolation order 1:
    phi.setorder(air, 1);
    
    // Specific weight of air under some circumstances:
    double gamma = 1.4;
    
    // Define the air density 'rho' and the Mach number:
    expression rho = pow(1 + (gamma-1)/2.0 * 0.7 * 0.7 * (1-grad(phi)*grad(phi)), 1.0/(gamma-1));
    expression machnumber = sqrt(grad(phi)*grad(phi)) / sqrt(1.0/(0.7*0.7) + 0.5*(gamma-1) * (1-grad(phi)*grad(phi)));
    
    // We suppose outside the air domain a uniform speed of 1 with direction left to right.
    // Since grad(phi) is the speed we have that grad(x) indeed gives us what we want.
    phi.setconstraint(gammaouter, x);
    
    // Define the potential flow formulation:
    formulation pf;
    
    // On the airfoil boundary the default Neumann condition dphi/dnormal = 0 
    // automatically ensures that there is no fluid entering the airfoil.
    // We thus do not need anything else than this:
    pf += integral(air, rho * grad(dof(phi)) * grad(tf(phi)) );
    
    // Start the nonlinear iteration with an all zero guess:
    vec sol(pf);
    
    double relres = 1;
    while (relres > 1e-7)
    {
        // Generate the formulation:
        pf.generate();
        // Get A and b to solve Ax = b:
        mat A = pf.A();
        vec b = pf.b();
        
        // Compute a relative residual:
        relres = (b - A*sol).norm() / b.norm();
        
        // Solve Ax = b:
        sol = solve(A, b);
        
        // Transfer the data from the solution vector to field 'phi' on the whole 'air' region:
        phi.setdata(air, sol);
        
        std::cout << "Current iteration has relative residual: " << relres << std::endl;
    }
    
    // Write the fluid speed (i.e. grad(phi)) and the Mach number:
    grad(phi).write(air, "flowspeed.pos", 1);
    machnumber.write(air, "machnumber.pos", 1);
    
    // Code validation line. Can be removed.
    std::cout << (machnumber.integrate(air, 3) < 62.4149 && machnumber.integrate(air, 3) > 62.4145);
}

