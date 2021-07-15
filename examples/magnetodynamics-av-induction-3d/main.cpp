// This example shows the resolution of a 3D magnetodynamic problem of magnetic induction using the
// so called a-v formulation. An aluminium tube is surrounded by 5 turns of a thick copper wire (1cm radius).
// Because of the alternating voltage (50 Hz) applied to the ends of the wire, eddy currents appear
// in the conducting aluminium tube. In this example the skin effect in the thick wire can also be observed.
//
// In order to remove the singularity (that comes from the magnetic equations) in the generated algebraic matrix
// a gauge condition is used. For that a spanning tree is created on the whole mesh, starting the growth on the
// region where the magnetic potential vector field 'a' will be constrained (here the domain boundary).
//
// This example was adapted from, and validated with an example developed for the GetDP 
// software (Patrick Dular and Christophe Geuzaine).


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'inductionheating.geo':
    int coil = 1, tube = 2, air = 3, coilskin = 4, tubeskin = 5, vin = 6, vout = 7, domainboundary = 8;
    
    // Load the mesh file. It can include curved elements.
    mesh mymesh("inductionheating.msh");
    
    // Define extra regions for convenience:
    int conductor = selectunion({coil,tube});
    int wholedomain = selectall();
    domainboundary = selectunion({domainboundary,vin,vout});
    
    parameter mu, sigma;
    // Define the magnetic permeability mu [H/m] everywhere (the materials here are not magnetic):
    mu|wholedomain = 4*getpi()*1e-7;
    // Conductivity of the copper coil and aluminium tube [S/m]:
    sigma|coil = 6e7; 
    sigma|tube = 3e7;
    
    // Set the working frequency to 50 Hz:
    setfundamentalfrequency(50);
    
    // Define a spanning tree to gauge the magnetic vector potential (otherwise the matrix to invert is singular).
    // Start growing the tree from the regions with constrained potential vector (here the domain boundary): 
    spanningtree spantree({domainboundary});
    
    // Use nodal shape functions 'h1' for the electric scalar potential 'v'.
    // Use edge shape functions 'hcurl' for the magnetic vector potential 'a'.
    // A spanning tree has to be provided to field 'a' for gauging.
    // Since the solution has a component in phase with the electric actuation
    // and a quadrature component we need 2 harmonics at 50Hz 
    // (harmonic 1 is DC, 2 is sine at 50Hz and 3 cosine at 50Hz).
    field a("hcurl", {2,3}, spantree), v("h1", {2,3});
    
    // Gauge the vector potential field on the whole volume:
    a.setgauge(wholedomain);
    
    // Select adapted interpolation orders for field a and v:
    a.setorder(wholedomain, 0);
    v.setorder(wholedomain, 1);
    
    // Put a magnetic wall (i.e. set field a to 0) all around the domain (no magnetic flux can cross it):
    a.setconstraint(domainboundary);
    // Also ground v on face 'vout':
    v.setconstraint(vout);
    // Set v to 1V on face 'vin' for the in-phase component and to 0 for the quadrature component:
    v.harmonic(2).setconstraint(vin, 1);
    v.harmonic(3).setconstraint(vin);
    
    // Define the weak magnetodynamic formulation:
    formulation magdyn;
    
    // The strong form of the magnetodynamic a-v formulation is 
    // 
    // curl( 1/mu * curl(a) ) + sigma * (dt(a) + grad(v)) = js, with b = curl(a) and e = -dt(a) - grad(v)
    //
    // Magnetic equation:
    magdyn += integral(wholedomain, 1/mu* curl(dof(a)) * curl(tf(a)) );
    magdyn += integral(conductor, sigma*dt(dof(a))*tf(a) + sigma* grad(dof(v))*tf(a) );
    // Electric equation:
    magdyn += integral(conductor, sigma*grad(dof(v))*grad(tf(v)) + sigma*dt(dof(a))*grad(tf(v)) );
    
    // Generate, solve and transfer the solution to fields a and v:
    magdyn.solve();
    
    // Write the magnetic induction field b = curl(a) [T], electric field e = -dt(a) - grad(v) [V/m] and current density j [A/m^2]:
    expression e = -dt(a) - grad(v);
    
    curl(a).write(wholedomain, "b.pos", 1);
    e.write(conductor, "e.pos", 1);
    (sigma*e).write(conductor, "j.pos", 1);
    v.write(conductor, "v.pos", 1);
    
    // Code validation line. Can be removed:
    double minj = norm(sigma*(-2*getpi()*50*a.harmonic(2) - grad(v.harmonic(3)))).min(tube,4)[0];
    std::cout << (minj < 216432 && minj > 216430);
}

