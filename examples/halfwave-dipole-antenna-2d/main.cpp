// This code simulates the electromagnetic wave radiation of a half-wave dipole antenna excited at 1 GHz.
// The simulation is performed in 2D (i.e. the fields are constant in the z direction). The goal of this code 
// is to provide a lightweight example that can be used with minimal modifications to simulate actual 3D antennas. 
//
// The radiation pattern of this dipole antenna can be obtained by calling 'interpolate' on the time-averaged power 
// expression for the points of a circle near the domain boundary and centered at the origin. The interpolated
// values can then be saved to .csv using 'writevector' and 'polarplot' could be used in Matlab for visualization. 


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    wallclock clk;
    
    // The domain regions as defined in 'dipoleantenna.geo':
    int air = 1, conductor = 2, feed = 3, boundary = 4;
    
    mesh mymesh("dipoleantenna.msh");
    
    // Define the whole non-conducting region for convenience:
    int wholedomain = regionunion({air, feed});
    
    // Edge shape functions 'hcurl' for the electric field E.
    // Because of the propagating EM waves the electric field has an 
    // in-phase and quadrature component (thus harmonics 2 and 3 are used):
    // E = E2 * sin(2*pi*1GHz*t) + E3 * cos(2*pi*1GHz*t).
    field E("hcurl", {2,3});
    field Es = E.harmonic(2), Ec = E.harmonic(3);
    
    // Use interpolation order 2 on the non-conducting region:
    E.setorder(wholedomain, 2);
    E.setorder(conductor, 1);
    
    // Define the speed of light [m/s] and the working frequency [Hz].
    // The total half-wave dipole antenna length is 15 cm. This corresponds to a 1GHz frequency:
    double c = 3e8, freq = 1e9;
    
    setfundamentalfrequency(freq);    
    
    // The dipole is a perfect conductor. We thus force E to 0 on it:
    E.setconstraint(conductor);
    // We force an electric field of 1 V/m in the y direction on the feed port 
    // for the sin component of E only (the cos component is forced to zero):
    Es.setconstraint(feed, array3x1(0,1,0));
    Ec.setconstraint(feed);
    
    formulation maxwell;
    
    // This is the weak formulation in time for electromagnetic waves:
    maxwell += integral(wholedomain, -curl(dof(E))*curl(tf(E)) - 1/(c*c)*dtdt(dof(E))*tf(E));
    
    // The OUTWARD pointing normal is required for the wave radiation condition.
    // Depending on the 'boundary' lines orientation the normal can be flipped.
    // To confirm the normal direction it is written to disk below:
    expression n = -normal(boundary);
    n.write(boundary, "normal.pos");
    
    // Silver-Mueller radiation condition to force outgoing electromagnetic waves:
    maxwell += integral(boundary, -1/c * crossproduct(crossproduct(n, dt(dof(E))), n) * tf(E)); 
    
    // Generate the algebraic matrices A and right handside b of Ax = b:
    maxwell.generate();
    
    // Get the solution vector x of Ax = b:
    vec solE = solve(maxwell.A(), maxwell.b());
    
    // Transfer the data in the solution vector to field E:
    E.setdata(wholedomain, solE);    
    
    clk.print("Resolution time (before post-processing):");
    
    // Save the electric field E and magnetic field H with an order 2 interpolation.
    // Since from Maxwell -mu0*dt(H) = curl(E) we have that H = 1/(mu0*w^2) * curl(dt(E)).
    double mu0 = 4*getpi()*1e-7;
    
    expression H = 1/(mu0*pow(2*getpi()*freq, 2)) * curl(dt(E));
    
    H.write(wholedomain, "H.pos", 2);
    E.write(wholedomain, "E.pos", 2);
    // Write the electric field at 50 timesteps of a period for a time visualization: 
    E.write(wholedomain, "E.pos", 2, 50);
    
    // Write the Poynting vector E x H.
    expression S = crossproduct(E, H);
    // This involves a product that makes new harmonics appear at 2x1GHz (harmonics 4 and 5) plus a constant component (harmonic 1).
    S.write(wholedomain, "S.pos", 2);
    // Also save in time at 50 timesteps of a period:
    S.write(wholedomain, "S.pos", 2, 50);
    
    // Code validation line. Can be removed.
    std::cout << (norm(crossproduct(Es,Ec)).max(wholedomain, 5)[0] < 0.00218416 && norm(crossproduct(Es,Ec)).max(wholedomain, 5)[0] > 0.00218413);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

