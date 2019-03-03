// This code simulates the static deflection and small-signal vibration of a collapsed 3D CMUT while supposing axisymmetry.
// A CMUT is a microscale ultrasonic transducer with an electrostatic actuation. In collapse mode the CMUT membrane is 
// pulled-in, i.e. touches the bottom of the cavity.
//
// The CMUT geometry and mesh in this example are created using sparselizard.
// The CMUT is made up of (from bottom layer to top):
//
// - a substrate (not in this example)
// - a ground electrode
// - an insulator layer (SiO2)
// - a vacuum cavity
// - a membrane layer (Silicon)
// - an electrode


#include "sparselizardbase.h"


using namespace mathop;

// Arguments are:
// 
// CMUT radius, thickness of the membrane, depth of the cavity, insulator thickness and x-length of the pillar supporting the membrane.
//
mesh createmesh(double r, double thmem, double thcav, double thins, double lpillar);

void sparselizard(void)
{	    
    // Define the CMUT geometric dimensions [m]:
    double r = 20e-6, thmem = 0.5e-6, thcav = 0.5e-6, thins = 0.3e-6, lpillar = 10e-6;
    
    // Axisymmetric assumption:
    setaxisymmetry();
    
    // The domain regions as defined in 'createmesh':
    int membrane = 1, pillar = 2, cavity = 3, ground = 4, electrode = 5, insulator = 6;
    
    // Create the geometry and the mesh:    
    mesh mymesh = createmesh(r, thmem, thcav, thins, lpillar);
    
    // Write the mesh for display:
    mymesh.write("cmutaxisym.msh");
    
    // Define additional regions:
    int contact = regionintersection({membrane, cavity});
    int solid = regionunion({membrane, pillar, insulator});
    int wholedomain = regionunion({solid,cavity});
    
    // Nodal shape functions 'h1' for v (the electric potential) and 
    // u (membrane displacement). Three components are used for u.
    // Field umesh is used to smoothly deform the mesh in the cavity.
    //
    field u("h1xyz"), v("h1"), umesh("h1xyz"), nodalforcebalance("h1xyz");
    
    // Clamp the insulator:
    u.setconstraint(insulator);
    
    v.setconstraint(electrode, 120);	
    v.setconstraint(ground, 0);	
    
    // Set a conditional constraint on compy(u), the y component of the mechanical deflection.
    // The constraint is active for all nodal degrees of freedom at which 'condexpr' is positive or zero. 
    // 
    // There is mechanical contact if either the deflection is larger than a given threshold 
    // (chosen as 95% of the gap size times 1+1e-6, not 100% to avoid the need of remeshing the cavity)
    // or if at the same time:
    // 
    // - compy(u) is greater than a given deflection (95% of the gap size)
    // - the y direction force balance on the node is negative, i.e. the node is pulled downwards
    //
    // Define conditional expressions for the conditions above. The value of a conditional 
    // expression is the second argument expression at all positions in space where the first 
    // argument is positive or zero. The value is the third argument expression otherwise.
    //
    expression nodeforcebalancecondition = expression(compy(nodalforcebalance), -1, 1);
    expression contactcondition(-compy(u)-thcav*0.95, nodeforcebalancecondition, -1);
    expression condexpr(-compy(u)-(thcav+1e-6*thcav)*0.95, 1, contactcondition);
    // Set the conditional constraint on the y component of the deflection field u:
    u.compy().setconditionalconstraint(contact, condexpr, -0.95*(thcav+1e-8*thcav));
    
    // Young's modulus [Pa], Poisson's ratio [], the density [kg/m^3] and the electric permittivity [F/m]:
    parameter E, nu, rho, epsilon;
    
    E|solid = 150e9; E|insulator = 66e9;
    nu|solid = 0.25; nu|insulator = 0.17;
    rho|solid = 2320; rho|insulator = 2200;
    
    epsilon|cavity = 8.854e-12;
    epsilon|solid = 11.7*8.854e-12;
    epsilon|insulator = 3.9*8.854e-12;

    // An electrostatic formulation is used for the electric problem.
    // An elasticity formulation is used for the mechanical problem.
    // Geometrical nonlinearity is taken into account.
    formulation electrostatics, elasticity;
	
    // Weak electrostatic formulation, computed on the mesh deformed by field umesh:
    electrostatics += integral(wholedomain, umesh, epsilon*grad(dof(v))*grad(tf(v)));
    
    // Weak elasticity formulation with geometrical nonlinearity 
    // (refer to the documentation to add prestress with the last argument).
    elasticity += integral(solid, predefinedelasticity(dof(u), tf(u), u, E, nu, 0.0));

    // Electrostatic forces, computed on the elements of the whole electric domain
    // but with mechanical deflection test functions tf(u) only on solid.
    // 
    // The electrostatic forces often appear in MEMS simulations and are thus predefined.
    // The inputs are the gradient of the test function of u defined on the mechanical domain,
    // the gradient of the previously computed electric potential field and the electric permittivity.
    //
    // The electrostatic forces are computed on the mesh deformed by field umesh.
    elasticity += integral(wholedomain, umesh, predefinedelectrostaticforce(tf(u,solid), grad(v), epsilon));
    
    
    // Solve the Laplace equation in the cavity to smoothly deform the mesh.
    // umesh is forced to field u on region solid:
    umesh.setconstraint(solid, u);

    formulation laplacian;
    laplacian += integral(wholedomain, grad(dof(compx(umesh)))*grad(tf(compx(umesh))) + grad(dof(compy(umesh)))*grad(tf(compy(umesh))) + grad(dof(compz(umesh)))*grad(tf(compz(umesh))));
    
    
    // NONLINEAR ITERATION TO GET THE STATIC DEFLECTION:
  	
    // Start with an all-zero solution vector for the elasticity formulation:
    vec solu(elasticity);
    
    double relresnorm = 1; int iter = 0; 
    while (relresnorm > 1e-8)
    {
        electrostatics.generate();
        vec solv = solve(electrostatics.A(), electrostatics.b());
        // Transfer the data from the solution vector to the v field:
        v.setdata(wholedomain, solv);
        // Write the electric field computed and saved on the geometry deformed by umesh.
        (-grad(v)).write(wholedomain, umesh, "E.pos");
        v.write(wholedomain, umesh, "v.pos");
        
        // Use the now known electric potential v to compute the membrane deflection:
        elasticity.generate();
        
        // Calculate the force balance on every node candidate for the conditional constraint. 
        // This requires the residual b-Ax (of the problem Ax = b) without the conditional constraints:
        elasticity.skipconditionalconstraints();
        // Argument 'true' is added so that the generated data is not cleared from the 'elasticity' object.
        vec fbvec = -(elasticity.b(true) - elasticity.A(true) * solu);
        // The conditional constraints are activated again:
        elasticity.skipconditionalconstraints(false);
        // The data in the 'fbvec' vector (containing the nodal force balance) is transferred to the
        // 'nodalforcebalance' field. '|u' is needed because 'nodalforcebalance' is not a 
        // field defined in the 'elasticity' formulation (field u is defined).
        nodalforcebalance.setdata(solid, fbvec|u);
        
        // Get the vector b and matrix A of problem Ax = b. Here they include the conditional constraints:
        vec b = elasticity.b();
        mat A = elasticity.A();
        // Compute the norm of the relative residual to check the convergence:
        relresnorm = (b-A*solu).norm()/b.norm();
        
        solu = solve(A,b);
        u.setdata(solid, solu);
        // Write the deflection u:
        u.write(solid, "u.pos"); 
        
        // Smooth the mesh in the cavity:
        laplacian.generate();
        vec solumesh = solve(laplacian.A(), laplacian.b());
        // Save the smoothed mesh in the cavity:
        umesh.setdata(wholedomain, solumesh);
        
        // Print the iteration number and relative residual norm:
        std::cout << "Rel. res. norm @it " << iter << " is " << relresnorm << ", max deflection is " << abs(compy(u)).max(contact, 5)[0]*1e9 << " nm" << std::endl;
        iter++;
    }
    
    
    // HARMONIC PERTURBATION AROUND THE STATIC DEFLECTION:
    
    // AC electric actuation frequency:
    setfundamentalfrequency(1e6);
    
    // Fields uh and vh have a constant deflection plus a harmonic deflection (harmonics 1 and 2):
    // uh = uh1 + uh2 * sin(2*pi*f0*t), vh = vh1 + vh2 * sin(2*pi*f0*t)
    field uh("h1xyz", {1,2}), vh("h1", {1,2});
    
    // Set the static deflection to the above solution:
    uh.harmonic(1).setdata(solid, solu|u);
    
    // Clamp the insulator:
    uh.setconstraint(insulator);
    
    // Set the same conditional constraint as above:
    uh.harmonic(1).compy().setconditionalconstraint(contact, condexpr, -0.95*(thcav+1e-8*thcav));
    // The vibration around the static deflection is 0 at the contact: 
    uh.harmonic(2).setconditionalconstraint(contact, condexpr, array3x1(0,0,0));
    
    // Set the DC voltage bias on the electrode:
    vh.harmonic(1).setconstraint(electrode, 100);	
    // Set a tiny AC voltage on vh2:
    vh.harmonic(2).setconstraint(electrode, 1);
    // Ground:	
    vh.setconstraint(ground);	
    
    
    // Redefine the electrostatic and elasticity formulations as above:
    formulation helectrostatics, helasticity;
    
    helectrostatics += integral(wholedomain, umesh, epsilon*grad(dof(vh))*grad(tf(vh)));
    helasticity += integral(solid, predefinedelasticity(dof(uh), tf(uh), u, E, nu, 0.0));
    // Here the electrostatic force must be computed using an FFT (on 5 timesteps)
    // because the force calculation involves nonlinear operations on multiharmonic fields:
    helasticity += integral(wholedomain, 5, umesh, predefinedelectrostaticforce(tf(uh,solid), grad(vh), epsilon));
    // The inertia term is added for the harmonic analysis:
    helasticity += integral(solid, -rho*dtdt(dof(uh))*tf(uh));
    
    // Generate and solve the electrostatic problem:
    helectrostatics.generate();
    vec solvh = solve(helectrostatics.A(), helectrostatics.b());
    vh.setdata(wholedomain, solvh);
    
    // Generate and solve the elasticity problem:
    helasticity.generate();
    vec soluh = solve(helasticity.A(), helasticity.b());
    uh.setdata(solid, soluh);
    uh.write(solid, umesh, "usmallsignal.pos");
    // Write the vibration at 50 timesteps of a period for a time visualization:
    (uh-uh.harmonic(1)).write(solid, umesh, "usmallsignal.pos", 1, 50);
    
    // Print the max AC deflection:
    double uacmax = abs(compy(uh.harmonic(2))).max(solid, 5)[0];
    std::cout << "Peak AC deflection: " << uacmax*1e9 << " nm" << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (uacmax < 0.977678e-9 && uacmax > 0.977675e-9);
}



// THE MESH BELOW IS FULLY STRUCTURED AND IS CREATED USING THE (BASIC) SPARSELIZARD GEOMETRY CREATION TOOL.
// THE ADVANTAGE OF IT IS THAT THE CODE ABOVE CAN BE CALLED FOR ANY CMUT DIMENSION WITHOUT NEEDING CALLS TO EXTERNAL MESHING SOFTWARE.

mesh createmesh(double r, double thmem, double thcav, double thins, double lpillar)
{
    // Give names to the physical region numbers:
    int membrane = 1, pillar = 2, cavity = 3, ground = 4, electrode = 5, insulator = 6;
    
    // Number of mesh layers:
    int nx = 40, nzcavity = 5, nzmembrane = 6, nzinsulator = 3;
    int nxpillar = nx*lpillar/r;
    
    // Insulator:
    double h = -thins;
    shape q01("quadrangle", insulator, {0,h,0, r,h,0, r,h+thins,0, 0,h+thins,0}, {nx, nzinsulator, nx, nzinsulator});
    shape q02("quadrangle", insulator, {r,h,0, r+lpillar,h,0, r+lpillar,h+thins,0, r,h+thins,0}, {nxpillar, nzinsulator, nxpillar, nzinsulator});
    
    // Cavity layer:
    h = 0.0;
    shape q11("quadrangle", cavity, {0,h,0, r,h,0, r,h+thcav,0, 0,h+thcav,0}, {nx, nzcavity, nx, nzcavity});
    shape q12("quadrangle", pillar, {r,h,0, r+lpillar,h,0, r+lpillar,h+thcav,0, r,h+thcav,0}, {nxpillar, nzcavity, nxpillar, nzcavity});
    
    // Ground line:
    shape groundline("union", ground, {q01.getsons()[0], q02.getsons()[0]});
    
    // Membrane layer:	
    h = h+thcav;
    shape q21("quadrangle", membrane, {0,h,0, r,h,0, r,h+thmem,0, 0,h+thmem,0}, {nx, nzmembrane, nx, nzmembrane});
    shape q22("quadrangle", membrane, {r,h,0, r+lpillar,h,0, r+lpillar,h+thmem,0, r,h+thmem,0}, {nxpillar, nzmembrane, nxpillar, nzmembrane});
    
    // Electrode line:
    shape electrodeline("union", electrode, {q21.getsons()[2], q22.getsons()[2]});
    
    // Provide to the mesh all shapes of interest:
    mesh mymesh({q01,q02,q11,q12,q21,q22,groundline,electrodeline});
    
    return mymesh;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

