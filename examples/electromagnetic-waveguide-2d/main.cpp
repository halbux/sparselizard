// This code simulates the steady-state electromagnetic waves in a cross-shaped 
// 2D waveguide made of a perfect conductor.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'waveguide3D.geo':
    int left = 1, skin = 2, wholedomain = 3;

    mesh mymesh("waveguide2D.msh");

    // Edge shape functions 'hcurl' for the electric field E.
    // Fields x and y are the x and y coordinate fields.
    field E("hcurl"), x("x"), y("y");

    // Use interpolation order 2 on the whole domain:
    E.setorder(wholedomain, 2);
    
    // The cutoff frequency for a 0.2 m width is freq = 0.75 GHz in theory. 
    // With this code and a fine enough mesh you will get the same value.
    double freq = 0.9e9, c = 3e8, pi = 3.14159, k = 2*pi*freq/c;
    
    // The waveguide is a perfect conductor. We thus force all
    // tangential components of E to 0 on the waveguide skin.
    E.setconstraint(skin);
    // We force an electric field in the y direction on region 'left'
    // that is 0 on the exterior of 'left' and one sine period inside.
    E.setconstraint(left, sin(y/0.1*pi)* array3x1(0,1,0));

    formulation maxwell;
    
    // This is the weak formulation for electromagnetic waves:
    maxwell += integral(wholedomain, -curl(dof(E))*curl(tf(E)) + k*k*dof(E)*tf(E));
    
    maxwell.generate();
    vec solE = solve(maxwell.A(), maxwell.b());
    
    E.setdata(wholedomain, solE);    
    // Save the electric field E and magnetic field H with an order 2 interpolation:
    curl(E).write(wholedomain, "H.pos", 2);
    E.write(wholedomain, "E.pos", 2);

    // Code validation line. Can be removed.
    std::cout << (solE.norm() < 0.6826 && solE.norm() > 0.6825);
}

