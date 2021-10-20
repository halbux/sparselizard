// This code computes on a 3D cube the electrostatic potential when 
// the left face is forced at 10 V and the right face is at 2 V.
// It is used to validate the AMR algorithm in 3D.


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
            vals[i] = std::pow(rnd, 28);
        }
    }

    return {dm};
}

int main(void)
{	
    slmpi::initialize();
    
    int rank = slmpi::getrank();

    // The domain regions as defined in 'cube.geo':
    int vol = 1, left = 2, right = 3, bnd = 4;
    
    std::string partname = allpartition("cube.msh");
    
    mesh mymesh;
    if (slmpi::isavailable())
        mymesh = mesh("gmsh:"+partname, bnd, 1);
    else
        mymesh = mesh("gmsh:"+partname);
    
    // Adaptive mesh refinement based on a random value on each element:
    expression randexpr(1, 1, randfunction, {0});
    mymesh.setadaptivity(randexpr, 0, 3);
    
    // Nodal shape functions 'h1' for the electric potential field:
    field v("h1");

    // Use interpolation order 2 on the whole domain:
    v.setorder(vol, 2);
    
    // Force 10 V on the left and 2 V on the right:
    v.setconstraint(left, 10);
    v.setconstraint(right, 2);
    
    // epsilon is the electric permittivity:
    double epsilon = 8.854e-12;
  
    formulation electrostatics;

    electrostatics += integral(vol, -epsilon*grad(dof(v))*grad(tf(v)));
    
    int numloops = 20;
    for (int i = 0; i < numloops; i++)
    {
        std::cout << "Iteration " << i << "/" << numloops << std::endl;
        electrostatics.allsolve(1e-15, 100, "lu", 0);
        
        alladapt(1);
     
        // Average slope:
        double slope = norm(grad(v)).integrate(vol, 5)/expression(1).integrate(vol, 5);
        double relerr = std::abs(slope-8.0)/8.0;
        if (relerr > 2e-14)
        {
            std::cout << "Error in AMR validation example - rel error was " << relerr << std::endl;
            std::cout << false;
            slmpi::finalize();
            return 0; 
        }
    }
    
    v.write(vol, "v"+std::to_string(1000+rank)+".vtu", 1);
    (-grad(v)).write(vol, "E"+std::to_string(1000+rank)+".vtu", 1);
    
    std::cout << true;
    slmpi::finalize();
    return 0; 
}

