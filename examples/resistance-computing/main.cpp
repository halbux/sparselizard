// This code shows the most accurate way to compute general resistance values.
// At the same time it demonstrates how to connect a lumped resistor to a
// FEM model of a resistor and how to feed the FEM model with a total current.
//
//      DC current
//       flow FEM
//       --------
//      |        |      R
// I -> |        | --/\/\/\/-- V0
//      |        |
//       --------


#include "sparselizard.h"


using namespace sl;

int main(void)
{
    // The domain regions as defined in 'quad.geo':
    int all = 1, left = 2, right = 3;
    
    mesh mymesh("quad.msh");
    
    // Nodal shape functions 'h1' for the electric potential field v:
    field v("h1"), x("x"), y("y");

    // Use interpolation order 2:
    v.setorder(all, 2);
    
    // Constraints on v are set with ports here
    
    // Electric conductivity [S/m] increasing with the y coordinate:
    expression sigma = 0.1 * (1+8*y);
    
    double R = 3; // Ohm
    
    // Define the V0 port not associated to the FEM problem: 
    port V0;
    
    // Associate a V/I (voltage/current) port pair to field v on the left
    // and right electrodes. Field v will be constant on each electrode.
    port Vl, Il, Vr, Ir;
    
    v.setport(left, Vl, Il);
    v.setport(right, Vr, Ir);
    //               |   |
    //     primal port   dual port
    //
    // The dual port holds the global Neumann term on the port region.
    // For an electrokinetic formulation this equals the total current.
  
    formulation electrokinetic;
    
    // Set a 1 A current flowing in through the left electrode:
    electrokinetic += Il - 1.0;
    // Link the right electrode to port V0:
    electrokinetic += Vr - (V0 - R*Ir);
    // Set V0 to 1V:
    electrokinetic += V0 - 1.0;

    // Define the weak formulation for the DC current flow:
    electrokinetic += integral(all, -sigma * grad(dof(v)) * grad(tf(v)));

    // Generate, solve and transfer the solution to field v and to the ports:
    electrokinetic.solve();
    
    double vl = Vl.getvalue(), vr = Vr.getvalue();
    double il = Il.getvalue(), ir = Ir.getvalue();
    
    // Compute the FEM resistance:
    double Rfem = (vl - vr)/il;
    
    std::cout << "Vl = " << vl << " V | Vr = " << vr << " V" << std::endl;
    std::cout << "Il = " << il << " A | Ir = " << ir << " A" << std::endl;
    std::cout << "FEM resistance is " << Rfem << " Ohm" << std::endl;
    
    // Write v and j for illustration:
    v.write(all, "v.pos", 2);
    (sigma * -grad(v)).write(all, "j.pos", 2);
    
    // Code validation line. Can be removed.
    std::cout << (Rfem < 2+1e-14 && Rfem > 2-1e-14 && vr < 4+1e-14 && vr > 4-1e-14);
}

