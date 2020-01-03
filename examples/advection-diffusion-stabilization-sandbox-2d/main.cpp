// This is a sandbox example to get a feeling with the stabilization methods
// proposed in sparselizard for advection-dominated advection-diffusion problems.
//
// The geometry is a square of side 'a' with 'n' mesh points in the x and y direction.


#include "sparselizardbase.h"


using namespace mathop;

mesh createmesh(double a, int n);

void sparselizard()
{	
    // Square surface and its four sides:
    int sur = 1, left = 2, right = 3, up = 4, down = 5;

    mesh mymesh = createmesh(1.0, 101);

    // Field c can be a chemical concentration:
    field c("h1"), x("x"), y("y");

    // Use interpolation order 2:
    c.setorder(sur, 2);

    // Steady supply of chemical concentration from the left side: 
    c.setconstraint(left, 1);

    // A 0.1 m/s fluid flow in the x direction is considered:
    expression v = array2x1(0.1,0);

    // Initial all zero concentration except in the square 
    // delimited by [0,0.25] for x and [0.4,0.6] for y:
    expression cinit = ifpositive(x-0.25, 0.0, ifpositive(abs(y-0.5)-0.1, 0.0, 1.0) );

    c.setvalue(sur, cinit);
    c.write(sur, "cinit.vtu", 2);

    // Tuning factor for the stabilization:
    double delta = 0.1;
    
    // End time for the simulation:
    double tend = 5.0;	
    
    // This will hold the solution vectors for every timestep:
    std::vector<vec> sol;


    std::cout << std::endl << "No stabilization:" << std::endl;
    formulation diffnostab;

    diffnostab += integral(sur, predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    // No normal gradient on horizontal sides:
    diffnostab += integral(up, compy(grad(dof(c)))*tf(c));
    diffnostab += integral(down, compy(grad(dof(c)))*tf(c));

    c.setvalue(sur, cinit);
    vec init(diffnostab);
    init.setdata(sur, c);
    genalpha genanostab(diffnostab, init, vec(diffnostab), vec(diffnostab), {true,true,true,true});
    genanostab.setparameter(0.75);
    sol = genanostab.runlinear(0, 0.1, tend, 2)[0];
    c.setdata(sur, sol[sol.size()-1]);
    c.write(sur, "cnostab.vtu", 2);


    std::cout << std::endl << "Isotropic:" << std::endl;
    formulation diffiso;

    diffiso += integral(sur,predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    // No normal gradient on horizontal sides:
    diffiso += integral(up,compy(grad(dof(c)))*tf(c));
    diffiso += integral(down,compy(grad(dof(c)))*tf(c));
    diffiso += integral(sur, predefinedstabilization("iso", delta, c, v, 0.0, 0.0));

    c.setvalue(sur, cinit);
    vec initiso(diffiso);
    initiso.setdata(sur, c);
    genalpha genaiso(diffiso, initiso, vec(diffiso), vec(diffiso), {true,true,true,true});
    genaiso.setparameter(0.75);
    sol = genaiso.runlinear(0, 0.1, tend, 2)[0];
    c.setdata(sur, sol[sol.size()-1]);
    c.write(sur, "ciso.vtu", 2);


    std::cout << std::endl << "Streamline anisotropic:" << std::endl;
    formulation diffaniso;

    diffaniso += integral(sur,predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    // No normal gradient on horizontal sides:
    diffaniso += integral(up,compy(grad(dof(c)))*tf(c));
    diffaniso += integral(down,compy(grad(dof(c)))*tf(c));
    diffaniso += integral(sur, predefinedstabilization("aniso", delta, c, v, 0.0, 0.0));

    c.setvalue(sur, cinit);
    vec initaniso(diffaniso);
    initaniso.setdata(sur, c);
    genalpha genaaniso(diffaniso, initaniso, vec(diffaniso), vec(diffaniso), {true,true,true,true});
    genaaniso.setparameter(0.75);
    sol = genaaniso.runlinear(0, 0.1, tend, 2)[0];
    c.setdata(sur, sol[sol.size()-1]);
    c.write(sur, "caniso.vtu", 2);



    // Define the strong form resual to be used by the stabilization methods below:
    expression dofres = v*grad(dof(c));
    expression res = v*grad(c);

    // Tuning parameter for streamline and for crosswind:
    double deltas = 0.1, deltac = 0.001;

    std::cout << std::endl << "Streamline and crosswind combined:" << std::endl;
    formulation diffcomb;

    diffcomb += integral(sur,predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    // No normal gradient on horizontal sides:
    diffcomb += integral(up,compy(grad(dof(c)))*tf(c));
    diffcomb += integral(down,compy(grad(dof(c)))*tf(c));
    diffcomb += integral(sur, predefinedstabilization("supg", deltas, c, v, 1e-6, dofres));
    diffcomb += integral(sur, predefinedstabilization("cws", deltac, c, v, 1e-6, res));
    c.setvalue(sur, cinit);
    vec initcomb(diffcomb);
    initcomb.setdata(sur, c);

    genalpha genacomb(diffcomb, initcomb, vec(diffcomb), vec(diffcomb), {true,true,true,true});
    genacomb.setparameter(0.75);
    sol = genacomb.runlinear(0, 0.1, tend, 2)[0];
    c.setdata(sur, sol[sol.size()-1]);
    c.write(sur, "ccombined.vtu", 2);
    

    // Code validation line. Can be removed.
    double cinterpolated = c.interpolate(sur, {0.76,0.399,0.0})[0];
    std::cout << (cinterpolated < 0.128046 && cinterpolated > 0.128043);
}


mesh createmesh(double a, int n)
{
    // Give names to the physical region numbers:
    int sur = 1, left = 2, right = 3, up = 4, down = 5;

    shape square("quadrangle", sur , {0,0,0, a,0,0, a,a,0, 0,a,0}, {n,n,n,n});

    shape line1 = square.getsons()[0];
    line1.setphysicalregion(down);
    shape line2 = square.getsons()[1];
    line2.setphysicalregion(right);
    shape line3 = square.getsons()[2];
    line3.setphysicalregion(up);
    shape line4 = square.getsons()[3];
    line4.setphysicalregion(left);

    // Provide to the mesh all shapes of interest:
    mesh mymesh({square,line1,line2,line3,line4});

    return mymesh;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

