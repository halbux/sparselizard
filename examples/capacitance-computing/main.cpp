// This code shows the most accurate way to compute general capacitance values.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'capacitor.geo':
    int dielectric = 1, air = 2, electrode = 3, ground = 4;
    
    mesh mymesh("capacitor.msh");
    
    // Select the entire domain:
    int all = selectall();
    
    // Nodal shape functions 'h1' for the electric potential field v:
    field v("h1");

    // Use interpolation order 3:
    v.setorder(all, 3);
    
    // Ground:
    v.setconstraint(ground);
    
    // Electric permittivity:
    parameter epsilon;
    
    epsilon|air = 8.854e-12;
    epsilon|dielectric = 3.9 * 8.854e-12;
    
    // Associate a V/Q (voltage/charge) port pair to field v on the electrode.
    // The electrode is supposed to be a perfect conductor with a constant v.
    port V, Q;
    v.setport(electrode, V, Q);
    //                   |  |
    //         primal port  dual port
    //
    // The dual port holds the global Neumann term on the port region.
    // For an electrostatic formulation this equals the electrode charge.
  
    formulation electrostatics;
    
    // Set a 0.1 nC charge per unit depth on the electrode:
    electrostatics += Q - 0.1e-9;

    electrostatics += integral(all, -epsilon * grad(dof(v)) * grad(tf(v)));

    // Generate, solve and transfer the solution to field v and to the ports:
    electrostatics.solve();
    
    // Compute the capacitance:
    double C = Q.getvalue() / V.getvalue();
    
    std::cout << "Capacitance is " << C << " F per unit depth" << std::endl;
    std::cout << "Electrode voltage is " << V.getvalue() << " V" << std::endl;
    
    // Write v and E for illustration:
    v.write(all, "v.pos", 2);
    (-grad(v)).write(all, "E.pos", 2);
    
    // Code validation line. Can be removed.
    std::cout << (C < 1.30635e-10 && C > 1.30633e-10);
}

