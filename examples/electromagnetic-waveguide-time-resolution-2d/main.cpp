// This code simulates the electromagnetic wave propagation in a cross-shaped 2D
// waveguide made of a perfect conductor. A time simulation is performed with
// initial all-zero conditions.


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
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
    double freq = 0.9e9, c = 3e8, pi = 3.14159;
    
    // The waveguide is a perfect conductor. We thus force all
    // tangential components of E to 0 on the waveguide skin.
    E.setconstraint(skin);
    // We force an electric field in the y direction on region 'left'
    // that is 0 on the exterior of 'left' and one sine period inside.
    // The electric field varies in time at frequency freq.
    E.setconstraint(left, sin(y/0.1*pi)* sin(2*pi*freq*t()) *array3x1(0,1,0));

    formulation maxwell;
    
    // This is the weak formulation for electromagnetic waves in time:
    maxwell += integral(wholedomain, -curl(dof(E))*curl(tf(E)) - 1/(c*c)*dtdt(dof(E))*tf(E));
    
    // Define the generalized alpha object to time-solve formulation 'maxwell' 
    // with initial all zero time derivative vectors 'vec(maxwell)'.
    // The general system to solve in time is M*dtdtx + C*dtx + K*x = b.
    //
    // The last argument is a vector 'isconstant' telling if the:
    //
    // - excitation vector b is constant in time for isconstant[0]
    // - matrix K is constant in time for isconstant[1]
    // - matrix C is constant in time for isconstant[2]
    // - matrix M is constant in time for isconstant[3]
    //
    // You can set isconstant[0] to 'true' if in b only the constraints
    // are time dependent but the other excitation sources are not. 
    //
    // Setting properly the 'isconstant' vector can give a dramatic speedup
    // since it may avoid reassembling or allow reusing the matrix factorization.
    //
    genalpha ga(maxwell, vec(maxwell), vec(maxwell), 1, {true, true, true, true});
    
    double timestep = 0.05*1.0/freq;
    int numsteps = 400;
    
    settime(0);
    for (int ts = 0; ts < numsteps; ts++)
    {
        // Compute one timestep:
        ga.next(timestep);
        
        // Write with an order 2 interpolation and with the name of your choice:
        E.write(wholedomain, "E"+std::to_string(ts+100)+".pos", 2); 
    }
    
    // Code validation line. Can be removed.
    std::cout << (ga.gettimederivative()[0].norm() < 3.0665e+09 && ga.gettimederivative()[0].norm() > 3.0661e+09); 
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

