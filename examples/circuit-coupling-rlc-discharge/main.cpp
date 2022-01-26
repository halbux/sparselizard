// This code shows how to connect a resistor simulated with FEM to an RLC circuit.
// A transient RLC discharge is computed. The time values are written to text file.
//
//
//                                            DC current
//         Vc0                                 flow FEM
//        100 V                                --------
//          /     L 100m       R 360     I    |        |
//      |--/    --/\/\/\/-- --/\/\/\/-- ->- V |  R200  | -- 
//    __|__                                   |        |   |
//    _____ C 0.1u                             --------    |
//      |                                                  |
//      |__________________________________________________|
//
//
// Solution is i(t) = e^(-a*t) * B * sin(wd*t)
//
// with
//
// a = R/(2*L), w0 = 1/sqrt(L*C), wd = sqrt(w0^2 - a^2), B = V0/(wd*L) and V0 = 100.


#include "sparselizard.h"


using namespace sl;

int main(void)
{
    // The domain regions as defined in 'quad.geo':
    int all = 1, left = 2, right = 3;
    
    mesh mymesh("quad.msh");
    
    // Nodal shape functions 'h1' for the electric potential field v:
    field v("h1");

    // Use interpolation order 2:
    v.setorder(all, 2);
    
    // Ground the right FEM side:
    v.setconstraint(right);
    
    // Electric conductivity [S/m]:
    expression sigma = 0.005;
    
    // Associate a V/I (voltage/current) port pair to field v on the electrode.
    // The electrode is supposed to be a perfect conductor with a constant v.
    port Vc, I, V;
    v.setport(left, V, I);
    //              |  |
    //    primal port  dual port
    //
    // The dual port holds the global Neumann term on the port region.
    // For an electrokinetic formulation this equals the total current.
    //
    // I is positive if it enters the FEM block.
  
    double R = 360, L = 100e-3, C = 0.1e-6;
  
    formulation elec;
    
    // Set the RLC circuit equations:
    elec += Vc - (V + R*I + L*dt(I));
    elec += dt(Vc) + 1/C*I;

    // Define the weak formulation for the conduction current flow:
    elec += integral(all, -sigma * grad(dof(v)) * grad(tf(v)));

    // Initial capacitor voltage is 100 V:
    Vc.setvalue(100);
    
    // Time resolution object (all matrices can be reused):
    genalpha ga(elec, vec(elec), vec(elec), 0, {true, true, true, true});

    // Lowest high-frequency damping:
    ga.setparameter(1);

    // Store the Vc and I values at all timesteps:
    std::vector<double> Vcvals = {}, Ivals = {}, tvals = {};

    settime(0);
    while (gettime() < 2000e-6)
    {
        Vcvals.push_back(Vc.getvalue());
        Ivals.push_back(I.getvalue());
        tvals.push_back(gettime());
        
        ga.next(5e-6);
    }
    
    std::cout << "FEM resistance is " << V.getvalue()/I.getvalue() << " Ohm" << std::endl;
    
    // Write the vectors to text file:
    writevector("Vc.txt", Vcvals);
    writevector("I.txt", Ivals);
    writevector("t.txt", tvals);
    
    // Code validation line. Can be removed.
    double maxI = Ivals[27];
    double minI = Ivals[92];
    std::cout << (std::abs(maxI - 6.8700e-02)/maxI < 2e-4 && std::abs(minI + 2.7478e-02)/std::abs(minI) < 6e-4);
}

