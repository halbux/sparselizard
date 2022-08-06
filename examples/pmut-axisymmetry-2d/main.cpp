// This code simulates the harmonic mechanical deflection of a 3D PMUT while supposing axisymmetry.
// A PMUT is a microscale ultrasonic transducer with a piezoelectric actuation. In this example
// the PMUT vibrates in air but it is straightforward to use another medium (e.g. water).
//
// The PMUT geometry and mesh in this example are created using sparselizard but GMSH could have been used instead.
// In any case the goal is not to illustrate the geometry and mesh creation but rather the PMUT axisymmetric simulation.
// The PMUT is made up of (from bottom layer to top):
//
// - a substrate (not in this example)
// - a vacuum cavity
// - a membrane layer
// - a bottom metal layer
// - a piezo layer (sandwiched between the two metal layers)
// - a top metal layer
// - a fluid
//
// Refer to the following paper for more details:
//
// "Characterization of polymer-based piezoelectric micromachined ultrasound transducers for short-range gesture 
// recognition applications", P. Gijsenbergh, A. Halbach et al.
//
// Please note that since this PMUT is being fabricated the materials have been modified for confidentiality reasons.
// Standard literature properties or arbitrary values have been considered for the material properties.
// It should thus be no surprise that the resonance frequencies and deflections do not match the paper values.


#include "sparselizard.h"


using namespace sl;

// Arguments are:
//
// PMUT radius, radius of the fluid region around (should be large enough), thickness of the top metal, piezo layer and bottom metal, 
// thickness of the membrane, thickness of the cavity, electrode coverage (in %), x-length of the pillar under the membrane and wavelength in the fluid.
//
mesh createmesh(double r, double rfluid, double thtopelec, double thpiezo, double thbotelec, double thmem, double thcav, double cov, double lpillar, double wavelength);

int main(void)
{
    wallclock clk;
    
    // Driving frequency [Hz]:
    double f0 = 165.31e3;
    
    // Acoustic propagation speed c [m/s] and a scaling factor for numerical conditionning:
    double c = 340, scaling = 1e10;

    // Define the PMUT geometric dimensions [m]:
    double r = 300e-6, rfluid = 20e-3, thtopelec = 100e-9, thpiezo = 500e-9, thbotelec = 100e-9, thmem = 15e-6, thcav = 35e-6, cov = 0.67, lpillar = 300e-6;

    // Axisymmetric assumption:
    setaxisymmetry();

    // The domain regions as defined in 'createmesh':
    int piezo = 1, membrane = 2, topelec = 3, botelec = 4, pillar = 5, fluid = 6, clamp = 7, electrode = 8, fluidboundary = 9;

    // Create the geometry and the mesh:
    mesh mymesh = createmesh(r, rfluid, thtopelec, thpiezo, thbotelec, thmem, thcav, cov, lpillar, c/f0);

    // Write the mesh for display:
    mymesh.write("pmutaxisym.msh");

    // Define additional regions:
    int ground = selectintersection({piezo, botelec}, 1);
    int solid = selectunion({pillar, membrane, botelec, piezo, topelec});
    int isotropicsolid = selectunion({pillar, membrane, botelec, topelec});
    int pmuttop = selectintersection({topelec,fluid}, 1);

    // Harmonic analysis. Set the fundamental frequency [Hz]:
    setfundamentalfrequency(f0);

    // Nodal shape functions 'h1' for v (the electric potential), p (acoustic
    // pressure) and u (membrane displacement). Three components are used for u.
    // Use harmonic 2 and 3, i.e. u(x,t) = Us(x)*sin(2pi*f0*t) + Uc(x)*cos(2pi*f0*t)
    // for u and p and harmonic 2, i.e. v(x,t) = V(x)*sin(2pi*f0*t) for v.
    // The y coordinate field is also defined since it might be of interest.
    //
    field u("h1xyz",{2,3}), p("h1", {2,3}), v("h1",{2}), y("y");

    // Use interpolation order 2 for p and u and 1 for v:
    u.setorder(solid, 2);
    p.setorder(fluid, 2);
    v.setorder(piezo, 1);

    // Clamp on the clamp region:
    u.setconstraint(clamp);

    v.harmonic(2).setconstraint(electrode, 1);
    v.setconstraint(ground, 0);

    // Young's modulus [Pa], Poisson's ratio [] (for isotropic materials) and the density [kg/m^3]:
    parameter E, nu, rho;

    E|pillar = 3.1e9;
    E|membrane = 3.1e9;
    E|botelec = 100e9;
    E|topelec = 100e9;

    nu|pillar = 0.34;
    nu|membrane = 0.34;
    nu|botelec = 0.3;
    nu|topelec = 0.3;

    rho|pillar = 1300;
    rho|membrane = 1300;
    rho|botelec = 1000;
    rho|piezo = 1780;
    rho|topelec = 1000;
    rho|fluid = 1.2;
    
    // Acoustic attenuation alpha [Neper/m] at the considered frequency:
    expression alpha = dbtoneper(5.4);
    
    // Diagonal relative permittivity matrix for the piezo:
    expression K(3,3,{7.4,7.4,8.0});
    K = K * 8.854e-12;

    // Coupling matrix [C/m^2] in Voigt notation for the piezo (6 rows, 3 columns).
    // The matrix is expressed in the usual coordinates, i.e. with z being the layer growth direction.
    expression C(6,3,{0,0,0.024, 0,0,0.024, 0,0,-0.027, 0,-0.015,0, -0.015,0,0, 0,0,0});

    // Elasticity matrix [Pa] in Voigt notation for the piezo.
    // The matrix is expressed in the usual coordinates, i.e. with z being the layer growth direction.
    expression H(6,6, {3.8e9, 1.9e9,3.8e9, 0.9e9,0.9e9,1.2e9, 0,0,0,7e8, 0,0,0,0,7e8, 0,0,0,0,0,9e8});

    // Take into account that here y is the growth direction and not z.
    // Refer to the documentation to make sure you understand how to use 'rotate'!
    K.rotate(-90,0,0); C.rotate(-90,0,0,"K","RT"); H.rotate(-90,0,0,"K","KT");
    

    formulation pmutmodel;

    // Standard isotropic elasticity. dof(u) is the unknown field, tf(u) the test function field.
    pmutmodel += integral(isotropicsolid, predefinedelasticity(dof(u), tf(u), E, nu) );
    // All inertia terms:
    pmutmodel += integral(solid, -rho*dtdt(dof(u))*tf(u) );

    // The classical weak formulation for piezoelectricity below can be found e.g. in paper:
    //
    // "A new fnite-element formulation for electromechanical boundary value problems"

    // Define the mechanical equations of the problem in the piezo.
    // strain(u) returns the strain vector [exx,eyy,ezz,gyz,gxz,gxy] of u.
    pmutmodel += integral(piezo, -( H*strain(dof(u)) )*strain(tf(u)) -( C*grad(dof(v)) )*strain(tf(u)) );
    // Define the electrical equations of the problem in the piezo:
    pmutmodel += integral(piezo, ( transpose(C)*strain(dof(u)) )*grad(tf(v)) - ( K*grad(dof(v)) )*grad(tf(v)) );

    // The wave equation is solved in the fluid:
    pmutmodel += integral(fluid, predefinedacousticwave(dof(p), tf(p), c, alpha));
    // A Sommerfeld condition is used on the fluid boundary to have outgoing waves:
    pmutmodel += integral(fluidboundary, predefinedacousticradiation(dof(p), tf(p), c, alpha));

    // Elastoacoustic coupling terms are predefined. They consist in two terms.
    // The first term is the pressure applied by the fluid normal to the PMUT top.
    // The second term comes from Newton's law: a membrane acceleration creates a
    // pressure gradient in the fluid.
    //
    // To have a good matrix conditionning the pressure source is divided by
    // the scaling factor and to compensate it multiplies the pressure force.
    // This leads to the correct membrane deflection but a pressure field divided by the scaling factor.
    pmutmodel += integral(pmuttop, predefinedacousticstructureinteraction(dof(p), tf(p), dof(u), tf(u), c, rho, array3x1(0,1,0), alpha, scaling));

    // Generate, solve and transfer the solution to the u, p and v fields:
    pmutmodel.solve();

    // Write the deflection, pressure and electric potential to file (with an order 2 interpolation for p and u).
    u.write(solid, "u.vtk", 2);
    (scaling*p).write(fluid, "p.vtk", 2);
    v.write(piezo, "v.vtk", 1);
    // Write p with 50 timesteps for illustration:
    // (scaling*p).write(fluid, "p.vtk", 2, 50);

    // Output the peak deflection:
    std::cout << "Peak in-phase deflection:   " << 1e9*abs(u.compy().harmonic(2)).max(solid, 6)[0] << " nm" << std::endl;
    std::cout << "Peak quadrature deflection: " << 1e9*abs(u.compy().harmonic(3)).max(solid, 6)[0] << " nm" << std::endl;
    // Output the pressure at 'rfluid' meters above the PMUT center:
    double pressureabove = sqrt( pow(scaling*p.harmonic(2),2) + pow(scaling*p.harmonic(3),2) ).interpolate(fluid, {0,rfluid,0})[0];
    double peakpressure = sqrt( pow(scaling*p.harmonic(2),2) + pow(scaling*p.harmonic(3),2) ).max(fluid, 5)[0];
    std::cout << "Pressure at " << 1000*rfluid << " mm above PMUT center: " << pressureabove << " Pa" << std::endl;
    std::cout << "Peak pressure is " << peakpressure << " Pa" << std::endl;

    clk.print("Total computation time:");

    // Code validation line. Can be removed.
    std::cout << (pressureabove < 0.278699 && pressureabove > 0.278697);
}

// THE MESH BELOW IS FULLY STRUCTURED AND IS CREATED USING THE (BASIC) SPARSELIZARD GEOMETRY CREATION TOOL.
// THE ADVANTAGE OF IT IS THAT THE CODE ABOVE CAN BE CALLED FOR ANY PMUT DIMENSION WITHOUT NEEDING CALLS TO EXTERNAL MESHING SOFTWARE.
// AS AN ALTERNATIVE, GMSH COULD HAVE BEEN USED TO EASILY DEFINE THE GEOMETRY AND CREATE A DELAUNAY MESH IN THE FLUID.

mesh createmesh(double r, double rfluid, double thtopelec, double thpiezo, double thbotelec, double thmem, double thcav, double cov, double lpillar, double wavelength)
{
    // Give names to the physical region numbers:
    int piezo = 1, membrane = 2, topelec = 3, botelec = 4, pillar = 5, fluid = 6, clamp = 7, electrode = 8, fluidboundary = 9;

    // Number of mesh layers:
    int nxpmut = 10, nzthick = 10, nzthin = 5;
    int nxpillar = nxpmut*lpillar/r;
    
    // Calculate the number of fluid elements to have about 10 elements per wavelength:
    int nair = std::ceil(10.0 * rfluid/wavelength);
    
    // Cavity layer:
    double h = -(thcav+thmem+thbotelec+thpiezo+thtopelec);
    shape q13("quadrangle", pillar, {r,h,0, r+lpillar,h,0, r+lpillar,h+thcav,0, r,h+thcav,0}, {nxpillar, nzthick, nxpillar, nzthick});

    shape clampline1 = q13.getsons()[0];
    clampline1.setphysicalregion(clamp);
    shape clampline2 = q13.getsons()[1];
    clampline2.setphysicalregion(clamp);

    // Membrane layer:
    h = h+thcav;
    shape q21("quadrangle", membrane, {0,h,0, r*cov,h,0, r*cov,h+thmem,0, 0,h+thmem,0}, {int(nxpmut*cov), nzthick, int(nxpmut*cov), nzthick});
    shape q22("quadrangle", membrane, {r*cov,h,0, r,h,0, r,h+thmem,0, r*cov,h+thmem,0}, {int(nxpmut*(1-cov)), nzthick, int(nxpmut*(1-cov)), nzthick});
    shape q23("quadrangle", membrane, {r,h,0, r+lpillar,h,0, r+lpillar,h+thmem,0, r,h+thmem,0}, {nxpillar, nzthick, nxpillar, nzthick});
    // Bottom electrode layer:
    h = h+thmem;
    shape q31("quadrangle", botelec, {0,h,0, r*cov,h,0, r*cov,h+thbotelec,0, 0,h+thbotelec,0}, {int(nxpmut*cov), nzthin, int(nxpmut*cov), nzthin});
    shape q32("quadrangle", botelec, {r*cov,h,0, r,h,0, r,h+thbotelec,0, r*cov,h+thbotelec,0}, {int(nxpmut*(1-cov)), nzthin, int(nxpmut*(1-cov)), nzthin});
    shape q33("quadrangle", botelec, {r,h,0, r+lpillar,h,0, r+lpillar,h+thbotelec,0, r,h+thbotelec,0}, {nxpillar, nzthin, nxpillar, nzthin});
    // Piezo layer:
    h = h+thbotelec;
    shape q41("quadrangle", piezo, {0,h,0, r*cov,h,0, r*cov,h+thpiezo,0, 0,h+thpiezo,0}, {int(nxpmut*cov), nzthin, int(nxpmut*cov), nzthin});
    shape q42("quadrangle", piezo, {r*cov,h,0, r,h,0, r,h+thpiezo,0, r*cov,h+thpiezo,0}, {int(nxpmut*(1-cov)), nzthin, int(nxpmut*(1-cov)), nzthin});
    shape q43("quadrangle", piezo, {r,h,0, r+lpillar,h,0, r+lpillar,h+thpiezo,0, r,h+thpiezo,0}, {nxpillar, nzthin, nxpillar, nzthin});
    // Top electrode layer:
    h = h+thpiezo;
    shape q51("quadrangle", topelec, {0,h,0, r*cov,h,0, r*cov,h+thtopelec,0, 0,h+thtopelec,0}, {int(nxpmut*cov), nzthin, int(nxpmut*cov), nzthin});
    shape q52("quadrangle", topelec, {r*cov,h,0, r,h,0, r,h+thtopelec,0, r*cov,h+thtopelec,0}, {int(nxpmut*(1-cov)), nzthin, int(nxpmut*(1-cov)), nzthin});
    shape q53("quadrangle", topelec, {r,h,0, r+lpillar,h,0, r+lpillar,h+thtopelec,0, r,h+thtopelec,0}, {nxpillar, nzthin, nxpillar, nzthin});

    shape electrodeline = q51.getsons()[0];
    electrodeline.setphysicalregion(electrode);

    // Fluid:
    shape q61("quadrangle", fluid, {0,0,0, r*cov,0,0, r*cov,r+lpillar,0, 0,r+lpillar,0}, {int(nxpmut*cov), int(nxpmut*(r+lpillar)/r), int(nxpmut*cov), int(nxpmut*(r+lpillar)/r)});
    shape q62("quadrangle", fluid, {r*cov,0,0, r,0,0, r,r+lpillar,0, r*cov,r+lpillar,0}, {int(nxpmut*(1-cov)), int(nxpmut*(r+lpillar)/r), int(nxpmut*(1-cov)), int(nxpmut*(r+lpillar)/r)});
    shape q63("quadrangle", fluid, {r,0,0, r+lpillar,0,0, r+lpillar,r+lpillar,0, r,r+lpillar,0}, {nxpillar, int(nxpmut*(r+lpillar)/r), nxpillar, int(nxpmut*(r+lpillar)/r)});

    shape l1 = q61.getsons()[2];
    shape l2("line", -1, {r*cov,r+lpillar,0, r*cov,sqrt(pow(rfluid,2)-pow(r*cov,2)),0}, nair);
    shape l3("arc", -1, {r*cov,sqrt(pow(rfluid,2)-pow(r*cov,2)),0, 0,rfluid,0, 0,0,0}, int(nxpmut*cov));
    shape l4("line", -1, {0,r+lpillar,0, 0,rfluid,0}, nair);
    shape q71("quadrangle", fluid, {l1,l2,l3,l4});

    shape l5 = q62.getsons()[2];
    shape l6("line", -1, {r,r+lpillar,0, r,sqrt(pow(rfluid,2)-pow(r,2)),0}, nair);
    shape l7("arc", -1, {r,sqrt(pow(rfluid,2)-pow(r,2)),0, r*cov,sqrt(pow(rfluid,2)-pow(r*cov,2)),0, 0,0,0}, int(nxpmut*(1-cov)));
    shape q72("quadrangle", fluid, {l5,l6,l7,l2});

    shape l8 = q63.getsons()[2];
    shape l9("line", -1, {r+lpillar,r+lpillar,0, rfluid*sqrt(2)/2,rfluid*sqrt(2)/2,0}, nair);
    shape l10("arc", -1, {rfluid*sqrt(2)/2,rfluid*sqrt(2)/2,0, r,sqrt(pow(rfluid,2)-pow(r,2)),0, 0,0,0}, nxpillar);
    shape q73("quadrangle", fluid, {l8,l9,l10,l6});

    shape l11("line", -1, {r+lpillar,0,0, rfluid,0,0}, nair);
    shape l12("arc", -1, {rfluid,0,0, rfluid*sqrt(2)/2,rfluid*sqrt(2)/2,0, 0,0,0}, int(nxpmut*(r+lpillar)/r));
    shape l13 = q63.getsons()[1];
    shape q74("quadrangle", fluid, {l11,l12,l9,l13});

    // Fluid boundary:
    shape fb1 = q71.getsons()[2];
    fb1.setphysicalregion(fluidboundary);
    shape fb2 = q72.getsons()[2];
    fb2.setphysicalregion(fluidboundary);
    shape fb3 = q73.getsons()[2];
    fb3.setphysicalregion(fluidboundary);
    shape fb4 = q74.getsons()[1];
    fb4.setphysicalregion(fluidboundary);

    // Provide to the mesh all shapes of interest:
    mesh mymesh({q13,q21,q22,q23,q31,q32,q33,q41,q42,q43,q51,q52,q53,q61,q62,q63,q71,q72,q73,q74, clampline1,clampline2, electrodeline, fb1,fb2,fb3,fb4});

    return mymesh;
}

