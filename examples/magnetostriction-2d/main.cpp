// This code shows how to compute the elastic deformation of a choke inductor subject
// to static Maxwell and magnetostriction stresses. The choke is supposed thin enough
// to use a plane stress mechanical approximation. This example is based on an example
// developed by M. Rossi and J. Le Besnerais for the GetDP software. Please refer to paper
// 'Vibration Reduction of Inductors under Magnetostrictive and Maxwell Forces Excitation'
// for the simulation details. Further magnetostriction model details can be found in
// 'Magnetostriction Measurement by Using Dual Heterodyne Laser Interferometers'.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'choke.geo':
    int gaps = 1, core = 2, windingpos = 3, windingneg = 4, air = 5, clamp = 7, boundary = 8;
    
    // Load the mesh, refine it by splitting once and define the domain boundary:
    mesh mymesh;
    mymesh.split(1);
    mymesh.selectskin(boundary);
    mymesh.load("choke.msh");

    int all = selectall();
    int solid  = selectunion({core, gaps});
    int nonmag = selectunion({air, gaps, windingpos, windingneg});

    // Winding current [A] times number of turns:
    double In = 200 * 100;
    // Calculate the current density:
    double jsource = In / expression(1).integrate(windingpos, 5);

    // Nodal shape functions 'h1' for the z component of the vector potential.
    field az("h1");
    az.setorder(all, 2);

    // Put a magnetic wall at the outer air boundary:
    az.setconstraint(boundary);
    
    // In 2D the vector potential only has a z component:
    expression a = array3x1(0,0,az);

    // Vacuum magnetic permeability [H/m]: 
    double mu0 = 4.0*getpi()*1e-7;

    // Define the permeability in all regions.
    parameter mu;

    mu|core = 1000.0*mu0;
    mu|nonmag = mu0;

    formulation magnetostatics;

    // The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
    magnetostatics += integral(all, 1/mu* curl(dof(a)) * curl(tf(a)) );
    // Add the current density source jsource [A/m2] in the z direction:
    magnetostatics += integral(windingpos, -array3x1(0,0,jsource) * tf(a));
    magnetostatics += integral(windingneg, -array3x1(0,0,-jsource) * tf(a));

    magnetostatics.solve();
    
    expression b = curl(a);

    expression B = norm(b);
    expression bx = compx(b);
    expression by = compy(b);

    az.write(all, "az.pos", 1);
    b.write(all, "b.pos", 1);
    
    
    // Compute the elastic deformation due to the Maxwell and magnetostrictive stresses.
    
    // Load the measured magnetostriction data:
    spline measureddata("magnetostriction.txt");
    
    // Define the expression giving lambda tangent as a function of the magnetic induction b.
    // This internally uses a natural cubic spline interpolation of the loaded data samples.
    expression lT(measureddata, norm(b));
    lT = lT * 1e-6; // data loaded is in um/m
    expression lN = -0.5 * lT;
    
    // Mechanical displacement field:
    field u("h1xy");
    u.setorder(all, 2);
    
    // E is Young's modulus. nu is Poisson's ratio:
    parameter E, nu;
    
    E|core = 220e9;
    E|gaps = 15e9;
    
    nu|core = 0.3;
    nu|gaps = 0.3;
    
    u.setconstraint(clamp);
    

    // Magnetostrictive strain (in the xy coordinates) in Voigt form (xx,yy,xy):
    expression ems = 1/(B*B) * array3x1(lT*bx*bx + lN*by*by, lT*by*by + lN*bx*bx, bx*by*(lT-lN));
    // Plane stress elasticity tensor:
    expression H = E/(1-nu*nu) * array3x3(1,nu,0, nu,1,0, 0,0,0.5*(1-nu));
    
    
    formulation elasticity;
    
    elasticity += integral(solid, predefinedelasticity(dof(u), tf(u), E, nu, "planestress"));
    // Magnetostatic Maxwell stresses:
    elasticity += integral(all, predefinedmagnetostaticforce(tf(u, solid), b/mu, mu) );
    // Magnetostrictive stresses:
    elasticity += integral(core, H * ems * strain(tf(u)) );
    
    elasticity.solve();
    
    double probeval = norm(u).interpolate(solid, {0, 0.33438, 0})[0];
    
    std::cout << "U max is " << norm(u).max(solid, 5)[0] << " and U probe is " << probeval << " m" << std::endl;
    
    u.write(solid, "u.pos", 1);
    
    // Code validation line. Can be removed.
    std::cout << (std::abs(probeval - 1.127e-07)/probeval < 1e-4);
}

