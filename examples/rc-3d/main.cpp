// This code computes the resistance and capacitance in a 3D geometry made of a 
// conducting trace connected to a circular-shaped parallel-plate air capacitor.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    wallclock clk;
    
    // The domain regions as defined in 'rc.geo':
    int electrode = 1, ground = 2, conductor = 3, dielectric = 4; 
    
    // Load the mesh:
    mesh mymesh("rc.msh");
    
    // Define the whole domain for convenience:
    int wholedomain = selectunion({conductor, dielectric});
    
    // Nodal shape functions 'h1' for the electric potential field.
    // Harmonics 2 and 3 are respectively the in-phase and quadrature 
    // component of the electric potential field.
    field v("h1", {2,3});
    
    // The working frequency is 1 MHz:
    double freq = 1e6;
    setfundamentalfrequency(freq);
    
    // Use interpolation order 1 on the whole domain:
    v.setorder(wholedomain, 1);

    // Force 1 V in-phase on the electrode and 0 V on the ground:
    double appliedvoltage = 1;
    v.harmonic(2).setconstraint(electrode, 1);
    v.harmonic(3).setconstraint(electrode, 0);
    v.setconstraint(ground, 0);
    
    // sigma is the electric conductivity of the conductor [S/m]:
    double sigma = 1e3;
    // epsilon is the electric permittivity of the air [F/m]:
    double epsilon = 8.854e-12;
    
    formulation electrokinetics;
    
    // Define the weak formulation for the conduction current in 
    // the conductor and displacement current in the dielectric.
    //
    // The strong form is:
    //
    // div(sigma * grad(v)) + div(epsilon * grad(dt(v))) = 0
    //	
    // with E = -grad(v) and J = sigma*E in the conductor. 
    //
    electrokinetics += integral(conductor, sigma*grad(dof(v))*grad(tf(v)));
    electrokinetics += integral(dielectric, epsilon*grad(dt(dof(v)))*grad(tf(v)));
    
    electrokinetics.generate();
    
    // Solve the algebraic problem Ax = b to get the solution vector x:
    vec solv = solve(electrokinetics.A(), electrokinetics.b());
    
    // Transfer the data from the solution vector to the v field:
    v.setdata(wholedomain, solv);
    // Write the electric potential:
    v.write(wholedomain, "v.pos", 1);
    
    // Compute the total current flowing through the electrode face
    // in an alternative (but less accurate) way to using ports.
    //
    // Since the computation involves a gradient that has to be 
    // calculated in the volume (and not on the electrode face) 
    // one can not simply call (normal(conductor)*J).integrate(electrode,4)
    // since with this a surface gradient will be calculated.
    // 'on()' is called to force the evaluation in the volume.
    
    double Iinphase = (-normal(conductor) * on(conductor, sigma*(-grad(v.harmonic(2)))) ).integrate(electrode, 4);
    double Iquadrature = (-normal(conductor) * on(conductor, sigma*(-grad(v.harmonic(3)))) ).integrate(electrode, 4);
    double normI = sqrt(Iinphase*Iinphase + Iquadrature*Iquadrature);
    // The voltage and currents are known, thus R and C are known as well:
    double R = appliedvoltage/pow(normI,2) * Iinphase;
    double ZC = appliedvoltage/pow(normI,2) * Iquadrature;
    double C = 1/(ZC*2*getpi()*freq);
    
    std::cout << std::endl << "I = " << Iinphase << " + " << Iquadrature << "j A" << std::endl;
    std::cout << std::endl << "R = " << R << " Ohm. C = " << C << " F (parallel plate formula gives " << epsilon*getpi()*(300e-6*300e-6)/30e-6 << ")" << std::endl;
    std::cout << "Refine the mesh for a better match." << std::endl << std::endl;
    
    clk.print("Total computation time: ");
    
    // Code validation line. Can be removed.
    std::cout << (C < 8.34113e-14 && C > 8.34111e-14);
}

