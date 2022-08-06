// This code simulates the steady-state vibration of a 2D CMUT (capacitive micromachined
// ultrasonic transducer). A CMUT is a micromembrane, actuated with an electrostatic 
// force, that vibrates in a fluid (water here). Due to the small dimensions at play
// the electric excitation frequency is in the megahertz range.
// 
// In this code the electric excitation is a sine wave so that the vibration is
// periodic. Due to the electroelastic nonlinearity however the membrane vibration
// and the pressure field includes new harmonics. These harmonics are computed 
// with the multiharmonic resolution method in which all calculations are brought 
// back to the undeformed configuration as detailed in the paper
// "Steady-State, Nonlinear Analysis of Large Arrays of Electrically Actuated 
// Micromembranes Vibrating in a Fluid" (A Halbach, C Geuzaine).


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'cmut2d.geo':
    int insulator = 1, pillars = 2, vacuumgap = 3, membrane = 4, fluid = 5, ground = 6, clamp = 6, electrode = 7, fluidboundary = 8;
    
    // The mesh can be curved!
    mesh mymesh("cmut2d.msh");
    
    // Define the region for the mechanical problem:
    int solid = selectunion({insulator, pillars, membrane});
    // Define the region for the electric problem:
    int electricdomain = selectunion({insulator, pillars, membrane, vacuumgap});
    // Define the membrane top line:
    int membranetop = selectintersection({membrane, fluid}, 1);
    
    // Nodal shape functions 'h1' for the electric potential field v, acoustic pressure
    // field p and membrane deflection u (u has components x and y).
    // umesh smoothly deforms the mesh in the vacuum gap for a given membrane deflection.
    //
    // Multiharmonic fields are used. This means that the time-dependency of a field h 
    // is approximated by a truncated Fourier series:
    //
    // h = h1 + h2*sin(2pi*f0*t) + h3*cos(2pi*f0*t) + h4*sin(2 * 2pi*f0*t) + h5*cos(2 * 2pi*f0*t) + ... 
    //
    // where f0 is the fundamental frequency, t is the time variable and hi are the Fourier coefficients.
    // The second argument in the field declaration lists all harmonics considered in the truncation.
    // For example the displacement field u is supposed to be well approximated by 
    // u = u1 + u4*sin(2 * 2pi*f0*t) + u5*cos(2 * 2pi*f0*t) while in the pressure field the constant is removed.
    //
    field u("h1xy",{1,4,5}), umesh("h1xy",{1,4,5}), v("h1",{2,3}), p("h1",{4,5});
    
    // Set the fundamental frequency f0 to 6 MHz:
    setfundamentalfrequency(6e6);
    
    // Since this is a multiharmonic resolution on a mesh deformed by a multiharmonic 
    // displacement field the calculations must be brought back on the undeformed 
    // configuration with a variable change [x y z]deformed = [x y z]undeformed + umesh. 
    //
    // This is the Jacobian matrix for the change of variable:
    //
    expression jac = array2x2(1+dx(compx(umesh)),dx(compy(umesh)),  dy(compx(umesh)),1+dy(compy(umesh)));
    expression invjac = inverse(jac);
    expression detjac = determinant(jac);
    // The expressions above might appear multiple times in a formulation 
    // and reusing them while generating a formulation will give a speedup:
    invjac.reuseit(); detjac.reuseit();

    // Use interpolation order:
    //
    // - 3 for u on the membrane
    // - 2 for u on the pillars
    // - 1 elsewhere
    //
    // - 2 for the pressure field p
    //
    // - 1 for the electric potential v

    u.setorder(membrane, 3);
    u.setorder(vacuumgap, 3);
    u.setorder(pillars, 2);
    u.setorder(insulator, 1);
    
    p.setorder(fluid, 2);
    
    v.setorder(electricdomain, 1);
    
    umesh.setorder(solid, 1);
    umesh.setorder(electricdomain, 1);
    
    // Clamp and ground all harmonics (i.e. 0 valued-Dirichlet conditions for u and v):
    u.setconstraint(clamp);
    v.setconstraint(ground);
    // Force the electric potential on the electrode to 250*sin(2pi*f0*t).
    // First set all v harmonics to 0 then set harmonic 2 to 250:
    v.setconstraint(electrode);
    v.harmonic(2).setconstraint(electrode, 250);
  
    // E is Young's modulus. nu is Poisson's ratio. rho is the volumic mass.
    // epsilon is the electric permittivity. 
    //
    // The membrane is polysilicon, the insulator is silicon dioxyde.
    
    parameter E, nu, rho, epsilon;
    
    E|insulator =		70e9;
    E|pillars =			150e9;
    E|membrane =		150e9;
    
    nu|insulator =		0.17;
    nu|pillars =		0.3;
    nu|membrane =		0.3;
    
    rho|insulator =		2200;
    rho|pillars =		2330;
    rho|membrane =		2330;
    rho|fluid = 		1000;
    
    epsilon|vacuumgap =	8.854e-12;
    epsilon|insulator =	3.9*8.854e-12;
    epsilon|pillars =	11.7*8.854e-12;
    epsilon|membrane =	11.7*8.854e-12;
    
    // c is the speed of sound in the water fluid. scaling is a scaling factor.
    double c = 1484, scaling = 1e10;
    
    // An electrostatic formulation is used for the electric problem.
    // An elastoacoustic formulation is used for the fluid-mechanic problem.
    formulation electrostatics, elastoacoustic;
	
    // Weak electrostatic formulation brought back to the undeformed configuration.
    // Since this is a nonlinear multiharmonic formulation we add an integer argument.
    // This means an FFT will be used to get the 20 first harmonics in the nonlinear term.
    electrostatics += integral(electricdomain, 20, epsilon*transpose(invjac*grad(dof(v)))*invjac*grad(tf(v)) * detjac );
    
    // The linear elasticity formulation is classical and thus predefined:
    elastoacoustic += integral(solid, predefinedelasticity(dof(u), tf(u), E, nu, "planestrain"));
    // Add the inertia terms:
    elastoacoustic += integral(solid, -rho*dtdt(dof(u))*tf(u));
    
    // Electrostatic forces, computed on the elements of the whole electric domain
    // but with mechanical deflection test functions tf(u) only on solid.
    // 
    // The electrostatic forces often appear in MEMS simulations and are thus predefined.
    // The inputs are the gradient of the test function of u defined on the mechanical domain,
    // the gradient of the previously computed electric potential field and the electric permittivity.
    //
    // The electrostatic forces must be computed on the mesh deformed by field umesh.
    // The calculations are brought back to the undeformed configuration.
    elastoacoustic += integral(electricdomain, 20, predefinedelectrostaticforce(grad(tf(u,solid))*transpose(invjac), invjac*grad(v), epsilon) * detjac );
    
    // The wave equation is solved in the fluid:
    elastoacoustic += integral(fluid, -grad(dof(p))*grad(tf(p)) -1/pow(c,2)*dtdt(dof(p))*tf(p));
    // A Sommerfeld condition is used on the fluid boundary to have outgoing waves:
    elastoacoustic += integral(fluidboundary, -1/c*dt(dof(p))*tf(p));
    
    // Elastoacoustic coupling terms.
    // The first term is the forces applied by the fluid on the membrane top.
    // The second term is Newton's law: a membrane acceleration creates a 
    // pressure gradient in the fluid.
    //
    // To have a good matrix conditionning the pressure source is divided by 
    // the scaling factor and to compensate it multiplies the pressure force.
    // This leads to the correct membrane deflection but a pressure field divided by the scaling factor.
    elastoacoustic += integral(membranetop, -dof(p)*tf(compy(u)) * scaling);
    elastoacoustic += integral(membranetop, rho*dtdt(dof(compy(u)))*tf(p) / scaling);
    
    
    // Solve the Laplace equation in the vacuum gap to smoothly deform the mesh.
    // umesh is forced to field u on region solid:
    umesh.setconstraint(solid, u);

    formulation laplacian;
    laplacian += integral(vacuumgap, grad(dof(compx(umesh)))*grad(tf(compx(umesh))) + grad(dof(compy(umesh)))*grad(tf(compy(umesh))) );
    
    
    // NONLINEAR ITERATION TO GET THE DEFLECTION:
  	
    // Start with an all-zero solution vector for the elastoacoustic formulation:
    vec solup(elastoacoustic);
    
    double relresnorm = 1; int iter = 0; 
    while (relresnorm > 1e-5)
    {
        electrostatics.generate();
        vec solv = solve(electrostatics.A(), electrostatics.b());
        // Transfer the data from the solution vector to the v field:
        v.setdata(electricdomain, solv);
        
        // Use the now known electric potential v to compute the membrane deflection:
        elastoacoustic.generate();
        
        vec b = elastoacoustic.b();
        mat A = elastoacoustic.A();
        // Compute the norm of the relative residual:
        relresnorm = (b-A*solup).norm()/b.norm();

        solup = solve(A,b);
        u.setdata(solid, solup);
        p.setdata(fluid, solup);
        
        // Smooth the mesh on the vacuum gap:
        laplacian.generate();
        vec solumesh = solve(laplacian.A(), laplacian.b());
        umesh.setdata(vacuumgap, solumesh);

        // Also save the u field on region solid to umesh.
        // This is done by selecting field u with |u on the solu vector.
        umesh.setdata(solid, solup|u);
        
        // Print the iteration number and relative residual norm:
        std::cout << "Relative residual norm @iteration " << iter << " is " << relresnorm << std::endl;
        iter++;
    }
    
    // Write the electric field with an order 1 interpolation at 50 timesteps.
    (-grad(v)).write(electricdomain, "E.pos", 1, 50);
    
    // Write the deflection u with an order 3 interpolation:
    u.write(solid, "u.pos", 3); 
    // Also write it at 50 timesteps for a time view.
    u.write(solid, "u.pos", 3, 50);
    // Do the same with p (remember to multiply it by the scaling factor to get the actual pressure):
    (scaling*p).write(fluid, "p.pos", 3, 50);
    
    // Code validation line. Can be removed.
    std::cout << (solup.norm() < 0.00025261 && solup.norm() > 0.00025260);
}

