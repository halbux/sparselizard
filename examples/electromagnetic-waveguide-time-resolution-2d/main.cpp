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
    double freq = 0.7e9, c = 3e8, pi = 3.14159;
    
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
    
    // Define the Newmark object to time-solve formulation 'maxwell' 
    // with initial all zero solution vectors 'vec(maxwell)'.
    // The general system to solve with Newmark is M*dtdtx + C*dtx + K*x = b.
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
    // since it may avoid reassembling or allow reusing the LU factorisation.
    //
    newmark nm(maxwell, vec(maxwell), vec(maxwell), {true, true, true, true});
    
    // Run the Newmark time resolution from the time in the first argument
    // to the time in the third argument by timesteps given as second argument.
    // The last argument is optional (default is 1). When set to an integer 'n'
    // only one every n timesteps is added to the output vector.
    std::vector<vec> solvec = nm.runlinear(0, 0.05*1.0/freq, 20*1.0/freq, 2);
    
    // Now save all data in the 'solvec' vector (which contains the solution at every nth timestep).
    for (int ts = 0; ts < solvec.size(); ts++)
    {
        std::cout << ts << " ";
        // Transfer the data to the electric field:
        E.setdata(wholedomain, solvec[ts]);
        // Write with an order 2 interpolation and with the name of your choice:
        E.write(wholedomain, "E"+std::to_string(ts+100)+".pos",2); 
    }
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}










