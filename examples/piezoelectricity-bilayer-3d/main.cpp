// This code simulates the harmonic mechanical deflection of a 3D polysilicon - PZT microbilayer when 
// the PZT (a piezoelectric material) is sandwiched between two electrically actuated electrodes.


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
    // The domain regions as defined in 'bilayer.geo':
    int pztlayer = 1, polysiliconlayer = 2, electrode = 3, ground = 4, clamp = 6, freeside = 7, elecbox = 8, mecabox = 9;
    
    mesh mymesh("bilayer.msh");
    
    int wholedomain = regionunion({pztlayer, polysiliconlayer});
    
    // Harmonic analysis. Set the fundamental frequency [Hz]:
    setfundamentalfrequency(1e4);
    
    // Nodal shape functions 'h1' for v (the electric potential) and u 
    // (the membrane displacement). Three components are used for u.
    // Use harmonic 2, i.e. u(x,t) = U(x)*sin(2pi*f0*t) for u and 
    // v(x,t) = V(x)*sin(2pi*f0*t) for v.
    //
    field u("h1xyz",{2}), v("h1",{2});
    
    // Use interpolation order 2 on 'wholedomain' for u and 1 for v in the PZT layer:
    u.setorder(wholedomain, 2);
    v.setorder(pztlayer, 1);
    
    // Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);
    // To force a displacement and see the resulting electric potential you could use:
    // u.compx().setconstraint(freeside, 1e-8);
    // u.setconstraint(freeside, array3x1(1e-8,1.5e-8,0.8e-8));
    
    v.harmonic(2).setconstraint(electrode, 10);	
    v.setconstraint(ground, 0);	
    
    
    // Young's modulus, Poisson's ratio and the density of polysilicon:
    double E = 169e9, nu = 0.22, rhopolysi = 2320;
    
    // Diagonal relative permittivity matrix for PZT:
    expression K(3,3,{1704,1704,1433});
    K = K * 8.854e-12;
    
    // Coupling matrix [C/m^2] for PZT (6 rows, 3 columns):
    expression C(6,3,{0,0,-6.62, 0,0,-6.62, 0,0,23.24, 0,17.03,0, 17.03,0,0, 0,0,0});
    
    // Anisotropic Hooke's matrix [Pa] for PZT. Ordering is [exx,eyy,ezz,gyz,gxz,gxy] (Voigt form).
    // Only the lower triangular part (top to bottom and left to right) is provided since it is symmetric.
    expression H(6,6, {1.27e11, 8.02e10,1.27e11, 8.46e10,8.46e10,1.17e11, 0,0,0,2.29e10, 0,0,0,0,2.29e10, 0,0,0,0,0,2.34e10});
    
    // PZT density [kg/m^3]:
    double rhopzt = 7500;
    
    
    formulation piezo;
    
    // Standard isotropic elasticity in the polysilicon (not piezoelectric):
    piezo += integral(polysiliconlayer, predefinedelasticity(dof(u), tf(u), E, nu) );
    // Inertia term:
    piezo += integral(polysiliconlayer, -rhopolysi*dtdt(dof(u))*tf(u) );
    
    // The classical weak formulation for piezoelectricity below can be found e.g. in paper:
    //
    // "A new fnite-element formulation for electromechanical boundary value problems"
    
    // Define the mechanical equations of the problem in the piezo.
    // strain(u) returns the strain vector [exx,eyy,ezz,gyz,gxz,gxy] of u.
    piezo += integral(pztlayer, -( H*strain(dof(u)) )*strain(tf(u)) - ( C*grad(dof(v)) )*strain(tf(u)) );
    // Inertia term for PZT:
    piezo += integral(pztlayer, -rhopzt*dtdt(dof(u))*tf(u) );
    // Define the electrical equations of the problem in the piezo:
    piezo += integral(pztlayer, ( transpose(C)*strain(dof(u)) )*grad(tf(v)) - ( K*grad(dof(v)) )*grad(tf(v)) );
    
    piezo.generate();
    
    vec sol = solve(piezo.A(), piezo.b());
    
    // Transfer the data from the solution vector to the v and u fields:
    u.setdata(wholedomain, sol);
    v.setdata(pztlayer, sol);
    
    // Write the deflection and electric potential on the volume boundary.
    u.write(mecabox, "u.pos", 2);
    v.write(elecbox, "v.pos", 1);
    // Harmonic fields can be saved in time. Use 50 time steps:
    u.write(mecabox, "utime.pos", 2, 50);
    v.write(elecbox, "vtime.pos", 1, 50);
    
    // Code validation line. Can be removed.
    std::cout << (compz(u.harmonic(2)).integrate(electrode, 5) < -1.3622e-15 && compz(u.harmonic(2)).integrate(electrode, 5) > -1.3623e-15);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

