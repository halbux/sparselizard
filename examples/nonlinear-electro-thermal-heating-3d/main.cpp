// This code simulates the heating that arises when DC current flows through a tungsten conductor 
// with a thin central part. The DC current flow is calculated with an electrokinetic formulation
// based on applied electric potentials at both ends and the temperature field is calculated 
// using the heat equation formulation, including the Joule heating due to the current flow.
// The dependency of the tungsten resistivity and its thermal conductivity on the temperature is taken into 
// account. Because of this the problem is nonlinear and requires an iterative resolution to converge.


#include "sparselizard.h"


using namespace sl;

mesh createmesh(void);

int main(void)
{	
    // Region number 
    // - 1 corresponds to the whole volume
    // - 2 to the left face (input, electrically actuated)
    // - 3 to the right face (output, grounded)
    //
    int volume = 1, input = 2, output = 3;
    
    // Create the geometry and the mesh:    
    mesh mymesh = createmesh();
    
    // Write to file for visualization:
    mymesh.write("meshed.msh");

    // Voltage [V] applied between the input and output faces:
    double appliedvoltage = 0.2;
    
    // Create the electric potential field v and temperature field T (both with nodal shape functions h1):
    field v("h1"), T("h1");
    
    // Use an order 2 interpolation for both fields:
    v.setorder(volume, 2);
    T.setorder(volume, 2);
    
    // Apply the requested voltage to the input and ground the output:
    v.setconstraint(input, appliedvoltage);
    v.setconstraint(output);
    // Assume the temperature is fixed to 293 K on the input and output faces:
    T.setconstraint(input, 293);
    T.setconstraint(output, 293);
    
    // The tungsten resistivity 'rho' [Ohm*m] is a function of the temperature (5.60e-8 at 293 K with a 0.0045 [K^-1] temperature coefficient): 
    parameter rho, k;
    rho|volume = 5.6e-8*(1+0.0045*(T-293));
    // Its thermal conductivity varies from 173 to 118 [W/(m*K)] between 293 and 1000 K: 
    k|volume = 173-(173-118)*(T-293)/(1000-293);
    
    // Expression for the electric field E [V/m] and current density j [A/m^2]:
    expression E = -grad(v);
    expression j = 1/rho * E;

    
    // Define the weak formulation for the static current flow.
    //
    // The strong form is:
    //
    // div(1/rho * grad(v)) = 0
    //	
    // with E = -grad(v)
    //	
    formulation electrokinetics;
    
    electrokinetics += integral(volume, grad(tf(v))*1/rho*grad(dof(v)));
    
    // Define the weak formulation for the heat equation.
    //
    // Can be found at https://en.wikiversity.org/wiki/Nonlinear_finite_elements/Weak_form_of_heat_equation
    //
    // Strong form is r*cp*dT/dt - div(k*grad(T)) = volumic heat sources [W/m^3]
    //
    // where r [kg/m^3] is the density and cp [J/(kg*K)] is the specific heat capacity.
    //
    formulation heatequation;
    
    // Here the time derivative is zero (static simulation):
    // heatequation += integral(volume, r*cp*dt(dof(T))*tf(T));
    heatequation += integral(volume, grad(tf(T))*k*grad(dof(T)));
    // Add the current heating source:
    heatequation += integral(volume, -j*E*tf(T));


    
    // Start with a uniform 293 K temperature everywhere:
    T.setvalue(volume, 293);
    
    // Initial all-zero solution vector for the heat equation:
    vec solheat(heatequation);
    
    double relres = 1;
    while (relres > 1e-7)
    {
        // Compute the static current everywhere:
        electrokinetics.generate();
        // Get A and b to solve Ax = b:
        vec solv = solve(electrokinetics.A(), electrokinetics.b());
        // Transfer the data from the solution vector to the v field to be used for the heat equation below:
        v.setdata(volume, solv);
        
        // Deduce the temperature field everywhere:
        heatequation.generate();
        // Get A and b to solve Ax = b:
        mat Aheat = heatequation.A();
        vec bheat = heatequation.b();
        
        // Compute a relative residual:
        relres = (bheat - Aheat*solheat).norm() / bheat.norm();
        
        // Solve Ax = b:
        solheat = solve(Aheat, bheat);
        
        // Transfer the data from the solution vector to the T field to be used for the electrokinetics above:
        T.setdata(volume, solheat);
        
        std::cout << "Current iteration has relative residual: " << relres << std::endl;
    }

    // Compute the total current flowing through the input face
    // in an alternative (but less accurate) way to using ports.
    //
    // Since the computation involves a gradient that has to be 
    // calculated in the volume (and not on the input face) 
    // one can not simply call (normal(volume)*j).integrate(input,4)
    // since with this a surface gradient will be calculated.
    // 'on()' is called to force the evaluation in the volume.
    double I = (-normal(volume)*on(volume, j)).integrate(input, 4);
    // Compute the electric resistance R between input and output:
    double R = appliedvoltage/I;
    
    std::cout << std::endl << "Resistance is " << R << " Ohm. Current is " << I << " A" << std::endl;
    
    // Compute the peak temperature on the whole volume:
    double peaktemperature = T.max(volume,2)[0];
    std::cout << "Peak temperature is " << peaktemperature << " K" << std::endl << std::endl;
    
    // Write v, j and T:
    v.write(volume, "v.pos", 2);
    j.write(volume, "j.pos", 2);
    T.write(volume, "T.pos", 2);
    
    // Code validation line. Can be removed.
    std::cout << (peaktemperature < 618.446 && peaktemperature > 618.444 && I < 1851.30 && I > 1851.28);
}

mesh createmesh(void)
{
    // Give names to the physical region numbers:
    int volume = 1, input = 2, output = 3;
    
    // Define the x, y and z coordinate fields:
    field x("x"), y("y"), z("z");
    
    double length = 0.2, thickness = 0.01, width = 0.05;
    
    // Define the 2D base to be extruded:
    shape leftquad("quadrangle", -1, {-length/2,-width/2,-thickness/2, -0.5*length/2,-width/2,-thickness/2, -0.5*length/2,width/2,-thickness/2, -length/2,width/2,-thickness/2}, {8,8,8,8});
    shape rightquad("quadrangle", -1, {0.5*length/2,-width/2,-thickness/2, length/2,-width/2,-thickness/2, length/2,width/2,-thickness/2, 0.5*length/2,width/2,-thickness/2}, {8,8,8,8});
    shape centralquad("quadrangle", -1, {-0.5*length/2,-width/2,-thickness/2, 0.5*length/2,-width/2,-thickness/2, 0.5*length/2,width/2,-thickness/2, -0.5*length/2,width/2,-thickness/2}, {16,8,16,8});
    
    // Extrude the 2D base:
    shape leftblock = leftquad.extrude(volume, thickness, 3);
    shape rightblock = rightquad.extrude(volume, thickness, 3);
    shape centralblock = centralquad.extrude(volume, thickness, 3);

    // Make the central block less large in the middle:
    centralblock.move(array3x1(0, 15*(abs(x)-0.5*length/2)*y, 0));
    // Make it also a little thinner (in the z direction) in the middle:
    centralblock.move(array3x1(0, 0, 15*(abs(x)-0.5*length/2)*z));
    
    // Get the input and output faces:
    shape inputface = leftblock.getsons()[4];
    shape outputface = rightblock.getsons()[2];
    
    // Assign the physical region numbers as detailed above:
    inputface.setphysicalregion(input);
    outputface.setphysicalregion(output);
    
    // Load to the mesh all shapes important for the finite element simulation:
    mesh mymesh({leftblock,rightblock,centralblock,inputface,outputface});
    
    // The mesh can be written at any time to have a feedback while creating the geometry!
    // mymesh.write("meshed.msh");

    return mymesh;
}

