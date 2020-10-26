// This example shows the time-simulation of a 2D axisymmetric high-temperature superconductor (HTS):
// a superconducting tube is used to shield its inside from an externally applied magnetic field 
// that increases linearly over time. For the simulation the so called a-v magnetodynamic formulation 
// is used. The v part can be disregarded in the configuration of this example.
//
// The tube height is 210 mm, its thickness is 1 mm and its inner radius is 10 mm. 
//
// This example is based on and validated with the experimental and simulation work done by K. Hogan in his thesis 
// "Passive magnetic shielding with bulk high temperature superconductors under a non-uniform magnetic field".
//
// Details on the equations for HTS can also be found in "Numerical simulation of the magnetization of 
// high-temperature superconductors: a 3D finite element method using a single time-step iteration" (G. Lousberg)


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{   
    // The physical region numbers as defined in the mesh file:
    int tube = 1, air = 2, domainskin = 3;
    
    setaxisymmetry();
    
    mesh mymesh("tubeinair.msh");
    
    // Define the whole domain for convenience:
    int wholedomain = regionunion({tube, air});
    
    // Magnetic permeability of the air:
    double mu = 4.0*getpi()*1e-7;
    
    // In this axisymmetric problem the magnetic vector potential can
    // be expressed as a vector with only a non-zero z component:
    field az("h1"), x("x");
    
    // Use interpolation order 1:
    az.setorder(wholedomain, 1);
    
    expression a = array3x1(0,0,az);
    
    // Put a magnetic wall (i.e. set field a to 0) all 
    // around the domain (no magnetic flux can cross it):
    az.setconstraint(domainskin);
    
    // The magnetic vector potential a is rewritten as the sum of a known
    // source term and an unknown part to be determined: atotal = a + asource.
    // We want to apply an external induction field b of 1 mT in 
    // the y direction that linearly increases over time:
    expression bsource = array3x1(0,1e-3*t(),0);
    // Since b = curl(a) (in cylindrical coordinates) a valid choice for a is:
    expression asource = array3x1(0,0,-1e-3*0.5*t()*x);
    // And dt(a)/dt is:
    expression dtasource = array3x1(0,0,-1e-3*0.5*x);
    
    
    // The strong form of the magnetodynamic a-v formulation is 
    // 
    // curl( 1/mu * curl(a) ) + sigma * (dt(a) + grad(v)) = js, with b = curl(a) and e = -dt(a) - grad(v)
    //
    // Here grad(v) is zero because of axisymmetry and y direction bsource and the electric field becomes:
    //
    expression e = -dt(a)-dtasource;
    
    // The conductivity of the high temperature superconductor is modeled using 
    // a power law relation between the electric field and the current density:
    //
    // J = Jc/Ec^(1/n) * norm(E)^((1-n)/n) * E
    //    |_______________________________|
    //                sigma(E)
    //
    // where the critical current density Jc [A/m^2] and critical electric field Ec [V/m] are supposed 
    // constant here and the J(E) relation is isotropic. An extra 1e-11 is added to never divide by 0.
    //
    double n = 30.0, Ec = 1e-4, Jc = 1e8;
    
    expression sigma = Jc/( pow(Ec,1.0/n) ) * pow( norm(e) + 1e-11, (1.0-n)/n );
    
    
    // Define the weak magnetodynamic formulation:
    formulation magdyn;
    
    magdyn += integral(wholedomain, 1/mu*( curl(dof(a)) + bsource ) * curl(tf(a)) );
    // Add an extra odd integration degree for convergence:
    magdyn += integral(tube, sigma*( dt(dof(a)) + dtasource )*tf(a), +1 );
    
    
    // Start the implicit Euler time resolution with the field values and a zero time derivative:
    impliciteuler eul(magdyn, vec(magdyn), 0);
    
    // Tolerance on the nonlinear iteration:
    eul.settolerance(1e-4);
    
    // Use a 1 second timestep and run until 150 sec (150 timesteps):
    int numsteps = 150;
    double timestep = 1.0;
    std::vector<double> bcenter(numsteps);
    
    settime(0);
    for (int i = 0; i < numsteps; i++)
    {
        // Use a fixed-point nonlinear iteration at every timestep (unlimited number of nl iterations with -1).
        int numnlits = eul.next(timestep, -1);
        
        // Write a and b = curl(a) to disk:
        norm(a).write(wholedomain, "norma" + std::to_string(i) + ".vtu"); 
        norm(curl(a) + bsource).write(wholedomain, "normbtotal" + std::to_string(i) + ".vtu"); 
        
        // Output the b induction field [T] at the tube center to assess the shielding effectiveness.
        // Interpolate at a x coordinate slightly away from 0 to avoid NaN issues:
        bcenter[i] = norm(curl(a) + bsource).interpolate(wholedomain, {1e-10,0,0})[0];
        std::cout << "b source " << timestep*(i+1) << " mT --> b tube center " << bcenter[i]*1e3 << " mT (" << numnlits << " NL its)" << std::endl;
    }
    // Combine all timesteps in a ParaView .pvd file for convenience:
    grouptimesteps("normbtotal.pvd", "normbtotal", 0, eul.gettimes());
    
    // Write to file the field at the tube center for all timesteps:
    writevector("bcenter.csv", bcenter);
    
    // Code validation line. Can be removed:
    std::cout << (bcenter[numsteps-1] < 0.03746 && bcenter[numsteps-1] > 0.03744);
}

int main(void)
{   
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

