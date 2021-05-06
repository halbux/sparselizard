// This code simulates the harmonic mechanical deflection of a 3D monocrystalline silicon - PZT microbilayer when 
// the PZT (a piezoelectric material) is sandwiched between two electrically actuated electrodes.
//
// The anisotropic behavior of both silicon and PZT is fully taken into account.
// The silicon and PZT crystals are rotated to align with a selected orientation.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'bilayer.geo':
    int pztlayer = 1, siliconlayer = 2, electrode = 3, ground = 4, clamp = 6, freeside = 7, elecbox = 8, mecabox = 9;
    
    mesh mymesh("bilayer.msh");
    
    int wholedomain = selectunion({pztlayer, siliconlayer});
    
    // Harmonic analysis. Set the fundamental frequency [Hz]:
    setfundamentalfrequency(1e4);
    
    // Nodal shape functions 'h1' for v (the electric potential) and u 
    // (the mechanical displacement). Three components are used for u.
    // Use harmonic 2, i.e. u(x,t) = U(x)*sin(2pi*f0*t) for u and 
    // v(x,t) = V(x)*sin(2pi*f0*t) for v.
    //
    field u("h1xyz",{2}), v("h1",{2});
    
    // Use interpolation order 2 on 'wholedomain' for u and 1 for v in the PZT layer:
    u.setorder(wholedomain, 2);
    v.setorder(pztlayer, 1);
    
    // Clamp on surface 'clamp' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);
    
    v.harmonic(2).setconstraint(electrode, 10);	
    v.setconstraint(ground, 0);	
    
    
    // Silicon and PZT density [kg/m^3]:
    double rhosi = 2330, rhopzt = 7500;
    
    // Diagonal relative permittivity matrix for PZT:
    expression K(3,3,{1704,1704,1433});
    K = K * 8.854e-12;
    
    // Coupling matrix [C/m^2] for PZT (6 rows, 3 columns):
    expression C(6,3,{0,0,-6.62, 0,0,-6.62, 0,0,23.24, 0,17.03,0, 17.03,0,0, 0,0,0});
    
    // Anisotropic elasticity matrix [Pa] for silicon and PZT. Ordering is [exx,eyy,ezz,gyz,gxz,gxy] (Voigt form).
    // Only the lower triangular part (top to bottom and left to right) is provided since it is symmetric.
    expression Hsi(6,6, {194.5e9, 35.7e9,194.5e9, 64.1e9,64.1e9,165.7e9, 0,0,0,79.6e9, 0,0,0,0,79.6e9, 0,0,0,0,0,50.9e9});
    expression Hpzt(6,6, {1.27e11, 8.02e10,1.27e11, 8.46e10,8.46e10,1.17e11, 0,0,0,2.29e10, 0,0,0,0,2.29e10, 0,0,0,0,0,2.34e10});
    
    // Rotate the silicon and PZT crystals first by -30 degrees around y then by 45 degrees around z.
    // Refer to the documentation to make sure you understand how to use 'rotate'!
    Hsi.rotate(0,-30,45,"K","KT"); Hpzt.rotate(0,-30,45,"K","KT"); C.rotate(0,-30,45,"K","RT"); K.rotate(0,-30,45);
    
    
    formulation piezo;
    
    // Orthotropic elasticity in the silicon (not piezoelectric):
    piezo += integral(siliconlayer, predefinedelasticity(dof(u), tf(u), Hsi) );
    // Inertia term:
    piezo += integral(siliconlayer, -rhosi*dtdt(dof(u))*tf(u) );
    
    // The classical weak formulation for piezoelectricity below can be found e.g. in paper:
    //
    // "A new fnite-element formulation for electromechanical boundary value problems"
    
    // Define the mechanical equations of the problem in the piezo.
    // strain(u) returns the strain vector [exx,eyy,ezz,gyz,gxz,gxy] of u.
    piezo += integral(pztlayer, -( Hpzt*strain(dof(u)) )*strain(tf(u)) - ( C*grad(dof(v)) )*strain(tf(u)) );
    // Inertia term for PZT:
    piezo += integral(pztlayer, -rhopzt*dtdt(dof(u))*tf(u) );
    // Define the electrical equations of the problem in the piezo:
    piezo += integral(pztlayer, ( transpose(C)*strain(dof(u)) )*grad(tf(v)) - ( K*grad(dof(v)) )*grad(tf(v)) );
    
    piezo.solve();
    
    // Display the peak displacement:
    double umax = norm(u.harmonic(2)).max(wholedomain, 5)[0];
    std::cout << "Peak deflection is " << umax << " m" << std::endl;
    
    // Write the deflection and electric potential on the volume boundary.
    u.write(mecabox, "u.vtk", 2);
    v.write(elecbox, "v.vtk", 1);
    // Harmonic fields can be saved in time. Use 50 time steps:
    u.write(mecabox, "utime.vtk", 2, 50);
    v.write(elecbox, "vtime.vtk", 1, 50);
    
    // Code validation line. Can be removed.
    std::cout << (umax < 1.7386e-07 && umax > 1.7384e-07);
}

