#include "mesh.h"
#include "field.h"
#include "expression.h"
#include "formulation.h"
#include "vec.h"
#include "petsc.h"
#include "wallclock.h"


using namespace mathop;

void sparselizard(void)
{	
    int vol = 1, sur = 2, top = 3;
    
    // The mesh can be curved!
    mesh mymesh("circle.msh");
    
    field u("h1xyz");

    u.setorder(vol, 3);
    
    u.setconstraint(sur);
  
    parameter E, nu;
    E|vol = 150e9; nu|vol = 0.3;
  
    formulation elasticity;

    elasticity += integral(vol, predefinedelasticity(u, E, nu));
    elasticity += integral(vol, array1x3(0,0,-10)*tf(u));

    elasticity.generate();

    vec solu = solve(elasticity.A(), elasticity.b());

    u.getdata(vol, solu);
    u.write(top, "u.pos", 3);
    
}

int main(void)
{	
    PetscInitialize(0,{},0,0);

    sparselizard();

    PetscFinalize();

    return 0;
}









