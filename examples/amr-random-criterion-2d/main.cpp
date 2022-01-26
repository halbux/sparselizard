// This code computes on a 2D quad the electrostatic potential when 
// the left line is forced at 10 V and the right line is at 2 V.
// It is used to validate the AMR algorithm in 2D.


#include "sparselizard.h"


using namespace sl;

std::vector<densemat> randfunction(std::vector<densemat> evaledexprs)
{
    int nrows = evaledexprs[0].countrows();
    int ncols = evaledexprs[0].countcolumns();
    
    densemat dm(nrows, ncols);
    
    double* vals = dm.getvalues();
    
    for (int i = 0; i < nrows*ncols; i++)
    {
        if (i == 0)
            vals[i] = 1; // Guarantee a 1 h-adapt criterion range
        else
        {
            double rnd = getrandom();
            // Get a good mix of each AMR depth:
            vals[i] = std::pow(rnd, 15);
        }
    }

    return {dm};
}

int main(void)
{	
    slmpi::initialize();
    
    int rank = slmpi::getrank();

    // The domain regions as defined in 'quad.geo':
    int sur = 1, left = 2, right = 3, bnd = 4;
    
    std::string partname = allpartition("quad.msh");
    
    mesh mymesh;
    if (slmpi::isavailable())
        mymesh = mesh("gmsh:"+partname, bnd, 1);
    else
        mymesh = mesh("gmsh:"+partname);
    
    // Adaptive mesh refinement based on a random value on each element:
    expression randexpr(1, 1, randfunction, {0});
    mymesh.setadaptivity(randexpr, 0, 4);
    
    // Nodal shape functions 'h1' for the electric potential field:
    field v("h1");

    // Use interpolation order 2 on the whole domain:
    v.setorder(sur, 2);
    
    // Force 10 V on the left and 2 V on the right:
    v.setconstraint(left, 10);
    v.setconstraint(right, 2);
    
    // epsilon is the electric permittivity:
    double epsilon = 8.854e-12;
  
    formulation electrostatics;

    electrostatics += integral(sur, -epsilon*grad(dof(v))*grad(tf(v)));
    
    int numloops = 100;
    for (int i = 0; i < numloops; i++)
    {
        std::cout << "Iteration " << i << "/" << numloops << std::endl;
        electrostatics.allsolve(1e-15, 100, "lu", 0);
        
        alladapt(0);
     
        // Average slope:
        double slope = norm(grad(v)).integrate(sur, 5)/expression(1).integrate(sur, 5);
        double relerr = std::abs(slope-8.0)/8.0;
        if (relerr > 1e-14)
        {
            std::cout << "Error in AMR validation example - rel error was " << relerr << std::endl;
            std::cout << false;
            slmpi::finalize();
            return 0; 
        }
    }
    
    v.write(sur, "v"+std::to_string(1000+rank)+".vtu", 1);
    (-grad(v)).write(sur, "E"+std::to_string(1000+rank)+".vtu", 1);
    
    std::cout << true;
    slmpi::finalize();
    return 0; 
}

