// This code computes on a 1D line the electrostatic potential when 
// the left point is forced at 10 V and the right point is at 2 V.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'line.geo':
    int line = 1, left = 2, right = 3;
    
    mesh mymesh("line.msh");
    
    // Nodal shape functions 'h1' for the electric potential field:
    field v("h1");

    // Use interpolation order 1 on the whole domain:
    v.setorder(line, 1);
    
    // Force 10 V on the left and 2 V on the right:
    v.setconstraint(left, 10);
    v.setconstraint(right, 2);
    
    // epsilon is the electric permittivity:
    double epsilon = 8.854e-12;
  
    formulation electrostatics;

    electrostatics += integral(line, -epsilon*grad(dof(v))*grad(tf(v)));

    electrostatics.generate();

    vec solv = solve(electrostatics.A(), electrostatics.b());

    // Transfer the data from the solution vector to the v field:
    v.setdata(line, solv);
    // Write v:
    v.write(line, "v.pos", 1);
    
    // Code validation line. Can be removed.
    std::cout << (solv.norm() < 21.5964 && solv.norm() > 21.5962);
}

