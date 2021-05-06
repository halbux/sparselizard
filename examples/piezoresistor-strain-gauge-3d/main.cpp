// This example shows how to set up a piezoresistive - elasticity simulation for a given silicon crystal orientation.
// As the silicon cantilever deflects the conduction current adapts due to a piezoresistive-induced change of the resistivity.
// The example also shows how to perform a field interpolation between two non-matching meshes.
//
// The piezoresistivity matrix was obtained from paper "Review: Semiconductor Piezoresistance for Microsystems", A. Barlian et al.
//
// Taking care of geometric nonlinearity (not needed here) and prestress can be easily done by calling the appropriate 
// 'predefinedelasticity' function, creating a nonlinear iteration to solve the mechanical problem and using the 
// 'greenlagrangestrain' function instead of the 'strain' function to correctly compute the strains (see documentation and examples).


#include "sparselizard.h"


using namespace sl;

// First arguments give the geometry size, last arguments give the number of nodes for meshing.
mesh createmesh(double length, double width, double thsi, double tracetilt, double ltrace, double wtrace, double rtrace, double thtrace, int numtraces, int nx, int ny, int nz, int nltrace, int nwtrace, int nthtrace);

int main(void)
{	
    int silicon = 1, trace = 2, electrode = 3, ground = 4, clamp = 5;
    
    // Create the two non-matching meshes:
    mesh mymesh = createmesh(200e-6, 100e-6, 10e-6, 30, 60e-6, 2e-6, 1e-6, 2e-6, 6, 20, 10, 5, 10, 4, 3);
    mymesh.write("straingauge.msh");
    
    // Nodal shape functions 'h1' for v (the electric potential) and u 
    // (the mechanical displacement). Three components are used for u. Field u is 
    // the displacement field in the silicon, field ugauge in the strain gauge.
    field u("h1xyz"), ugauge("h1xyz"), v("h1");
    
    // Use interpolation order 2 for u and for v:
    u.setorder(silicon, 2);
    u.setorder(trace, 2);
    ugauge.setorder(trace, 2);
    v.setorder(trace, 2);
    
    // Clamp on surface 'clamp' (i.e. 0 valued-Dirichlet conditions):
    u.setconstraint(clamp);
    
    v.setconstraint(electrode, 10);	
    v.setconstraint(ground, 0);	
    
    // Piezoresistivity 'Pi' [1/Pa] matrix aligned with the silicon crystal.
    // In this case it is symmetric and thus only the lower triangular part is provided:
    expression Pi(6,6,{6.6, -1.1,6.6, -1.1,-1.1,6.6, 0,0,0,138.1, 0,0,0,0,138.1, 0,0,0,0,0,138.1});
    Pi = Pi * 1e-11;
    // Multiply Pi by the unstressed resistivity to have an equation for the resistivity variation (and not the relative change).
    // The unstressed resistivity is 7.8e-2 Ohm*m. Pi has now units m^4/(s*A^2) or Ohm*m/Pa.
    double rhounstressed = 7.8e-2;
    Pi = Pi * rhounstressed;
    
    // Bring the Pi matrix to a 45 degrees crystal orientation.
    // Refer to the documentation to make sure you understand how to use 'rotate'!
    Pi.rotate(0,0,45,"K","K-1");
    
    // Anisotropic elasticity matrix [Pa] for silicon. Ordering is [exx,eyy,ezz,gyz,gxz,gxy] (Voigt form).
    // Only the lower triangular part (top to bottom and left to right) is provided since it is symmetric.
    expression H(6,6, {194.5e9, 35.7e9,194.5e9, 64.1e9,64.1e9,165.7e9, 0,0,0,79.6e9, 0,0,0,0,79.6e9, 0,0,0,0,0,50.9e9});
    
    // Resistivity change for the stressed silicon crystal (piezoelectric):
    expression deltarho = Pi * H*strain(ugauge);
    // Add the unstressed diagonal resistivity tensor in Voigt form:
    expression rhos = deltarho + expression(6,1,{rhounstressed,rhounstressed,rhounstressed,0,0,0});
    
    // Bring the resistivity from Voigt to 3x3 tensor form for inversion.
    // Only the lower triangular part is provided since the tensor is here symmetric.
    expression rho = expression(3,3, {rhos.at(0,0),rhos.at(5,0),rhos.at(1,0),rhos.at(4,0),rhos.at(3,0),rhos.at(2,0)});
    // The conductivity sigma [S/m] is the inverse of the resistivity:
    expression sigma = inverse(rho);
    
    // Mechanical formulation:
    formulation elasticity;
    
    // Orthotropic elasticity in the silicon:
    elasticity += integral(silicon, predefinedelasticity(dof(u), tf(u), H) );
    // Add a volume force [N/m^3] to deflect the cantilever:
    elasticity += integral(silicon, -1e11*compz(tf(u)) );
    
    // Define the weak formulation for the conduction current flow.
    //
    // The strong form is:
    //
    // div(1/rho * grad(v)) = 0
    //	
    // with E = -grad(v)
    //	
    formulation currentflow;
    
    currentflow += integral(trace, grad(tf(v))*(sigma*grad(dof(v))));
    
    
    // First solve the mechanic problem to get the stresses:
    elasticity.solve();
    
    // Transfer the deflection u to the ugauge field using mesh-to-mesh interpolation:
    ugauge.setvalue(trace, on(silicon, u));
    
    // Solve the DC current flow problem:
    currentflow.solve();
    
    // Expression for the electric field E [V/m] and current density j [A/m^2]:
    expression E = -grad(v);
    expression j = sigma * E;
    
    // Write the fields to file:
    u.write(silicon, "u.vtk", 2);
    ugauge.write(trace, "ugauge.vtk", 2);
    v.write(trace, ugauge, "v.vtk", 2);
    
    
    // Output the total current [A] flowing through the electrode.
    // The normal current density must be integrated on the whole electrode surface.
    // Since j requires calculating volume derivatives its evaluation must be performed on the 'trace' volume region.
    double I = (-normal(trace)*on(trace, j)).integrate(electrode, 4);
    std::cout << "Input current = " << I << " A at a peak deflection of " << norm(u).max(silicon,3)[0]*1e6 << " um" << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (I < 6.10153e-07 && I > 6.10150e-07);
}

mesh createmesh(double length, double width, double thsi, double tracetilt, double ltrace, double wtrace, double rtrace, double thtrace, int numtraces, int nx, int ny, int nz, int nltrace, int nwtrace, int nthtrace)
{
    int silicon = 1, trace = 2, electrode = 3, ground = 4, clamp = 5;

    // Create the cantilever:
    shape footprint("quadrangle", -1, {0,0,0, length,0,0, length,width,0, 0,width,0}, {nx,ny,nx,ny});
    shape si = footprint.extrude(silicon, thsi, nz);
    
    shape clmp = si.getsons()[4];
    clmp.setphysicalregion(clamp);
    
    // Create the strain gauge:
    shape q("quadrangle", -1, {0,0,0, ltrace,0,0, ltrace,wtrace,0, 0,wtrace,0}, {nltrace,nwtrace,nltrace,nwtrace});
    shape ain("arc", -1, {ltrace,wtrace,0, ltrace,wtrace+2.0*rtrace,0, ltrace,wtrace+rtrace,0}, nltrace);
    shape aout("arc", -1, {ltrace,0,0, ltrace,2.0*wtrace+2.0*rtrace,0, ltrace,wtrace+rtrace,0}, nltrace);
    shape ln("line", -1, {ltrace,wtrace+2.0*rtrace,0, ltrace,2.0*wtrace+2.0*rtrace,0}, nwtrace);
    shape qcurve("quadrangle", -1, {q.getsons()[1],aout,ln,ain});
    
    shape oneside("union", -1, {q,qcurve});
    shape otherside = oneside.duplicate();
    otherside.scale(-1,1,1);
    otherside.shift(ltrace, wtrace+2.0*rtrace,0);
    
    shape tracepiece("union", -1, {oneside, otherside});
    
    std::vector<shape> pieces(numtraces);
    for (int i = 0; i < numtraces; i++)
    {
        pieces[i] = tracepiece.duplicate();
        pieces[i].shift(0,2.0*(wtrace+2.0*rtrace)*i,0);
    }
    shape tracefootprint("union", -1, pieces);
    // Center:
    tracefootprint.shift(-ltrace/2,-numtraces/2*(4*rtrace+2*wtrace),0);

    shape trc = tracefootprint.extrude(trace, thtrace, nthtrace);

    trc.rotate(0,0,tracetilt);
    trc.shift(length/2, width/2, thsi-thtrace);

    // Get the electrode and ground faces:
    shape el = trc.getsons()[0].getsons()[0].getsons()[0].getsons()[4];
    el.setphysicalregion(electrode);
    shape gr = trc.getsons()[numtraces-1].getsons()[1].getsons()[1].getsons()[3];
    gr.setphysicalregion(ground);

    
    mesh mymesh({si,clmp,trc,el,gr});

    return mymesh;
}

