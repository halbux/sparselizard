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
    
    // Timestep and end time for the simulation:
    double ts = 0.1, tend = 4.99;


    std::cout << std::endl << "No stabilization:" << std::endl;
    formulation adnostab;

    adnostab += integral(sur, predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    
    c.setvalue(sur, cinit);
    genalpha genanostab(adnostab, vec(adnostab), vec(adnostab), 1, {true,true,true,true});
    genanostab.setparameter(0.75);
    settime(0);
    while (gettime() < tend)
        genanostab.next(ts);
    c.write(sur, "cnostab.vtu", 2);
        

    std::cout << std::endl << "Isotropic:" << std::endl;
    formulation adiso;

    adiso += integral(sur, predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    adiso += integral(sur, predefinedstabilization("iso", delta, c, v, 0.0, 0.0));

    c.setvalue(sur, cinit);
    genalpha genaiso(adiso, vec(adiso), vec(adiso), 1, {true,true,true,true});
    genaiso.setparameter(0.75);
    settime(0);
    while (gettime() < tend)
        genaiso.next(ts);
    c.write(sur, "ciso.vtu", 2);
    

    std::cout << std::endl << "Streamline anisotropic:" << std::endl;
    formulation adaniso;

    adaniso += integral(sur, predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    adaniso += integral(sur, predefinedstabilization("aniso", delta, c, v, 0.0, 0.0));

    c.setvalue(sur, cinit);
    genalpha genaaniso(adaniso, vec(adaniso), vec(adaniso), 1, {true,true,true,true});
    genaaniso.setparameter(0.75);
    settime(0);
    while (gettime() < tend)
        genaaniso.next(ts);
    c.write(sur, "caniso.vtu", 2);


    // Define the strong form residual to be used by the stabilization methods below:
    expression dofres = v*grad(dof(c));
    expression res = v*grad(c);

    // Tuning parameter for streamline and for crosswind:
    double deltas = 0.1, deltac = 0.001;

    std::cout << std::endl << "Streamline and crosswind combined:" << std::endl;
    formulation adcomb;

    adcomb += integral(sur, predefinedadvectiondiffusion(dof(c), tf(c), v, 1e-6, 1.0, 1.0, true));
    adcomb += integral(sur, predefinedstabilization("supg", deltas, c, v, 1e-6, dofres));
    adcomb += integral(sur, predefinedstabilization("cws", deltac, c, v, 1e-6, res));
    
    c.setvalue(sur, cinit);
    genalpha genacomb(adcomb, vec(adcomb), vec(adcomb), 1, {true,true,true,true});
    genacomb.setparameter(0.75);
    settime(0);
    while (gettime() < tend)
        genacomb.next(ts);
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

