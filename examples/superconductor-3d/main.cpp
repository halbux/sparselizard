// This example shows the time-simulation of a 3D high-temperature superconductor (HTS):
// a superconducting tube is used to shield its inside from an externally applied magnetic field 
// that increases linearly over time. For the simulation the so called a-v magnetodynamic formulation is used.
//
// The tube height is 210 mm, its thickness is 1 mm and its inner radius is 10 mm. 
//
// This example is based on and validated with the experimental and simulation work done by K. Hogan in his thesis 
// "Passive magnetic shielding with bulk high temperature superconductors under a non-uniform magnetic field".
//
// Details on the equations for HTS can also be found in "Numerical simulation of the magnetization of 
// high-temperature superconductors: a 3D finite element method using a single time-step iteration" (G. Lousberg)


#include "sparselizardbase.h"


using namespace mathop;

// Arguments are:
// 
// Tube thickness, height and inner radius, size of the air domain around as well as integers defining the mesh size.
// The last argument gives the mesh progression in the air around the tube (increase for a finer mesh at the tube).
//
mesh createmesh(double thtube, double htube, double rintube, double linf, int nhtube, int ncircumtube, int nthtube, int nair, int prog);

void sparselizard(void)
{	
    wallclock clk;
    
    // The domain regions as defined in 'createmesh':
    int wholedomain = 1, tube = 2, insidetube = 3, tubeskin = 4, domainskin = 5, ground = 6;
    
    // Tube thickness, height, inner radius and domain size [m]:
    double thtube = 1e-3, htube = 210e-3, rintube = 10e-3, linf = 100e-3;
    
    // Create the geometry and the mesh:   
    mesh mymesh = createmesh(thtube, htube, rintube, linf, 7, 4, 10, 7, 20);
    
    // Write the mesh for display:
    mymesh.write("tubeinair.msh");
    
    
    // Magnetic permeability of the air:
    double mu = 4*getpi()*1e-7;
    
    // Define a spanning tree to gauge the magnetic vector potential (otherwise the matrix to invert is singular).
    // Start growing the tree from the regions with constrained potential vector (here the domain boundary): 
    spanningtree spantree({domainskin});
    
    // Use nodal shape functions 'h1' for the electric scalar potential 'v'.
    // Use edge shape functions 'hcurl' for the magnetic vector potential 'a'.
    // A spanning tree has to be provided to field 'a' for gauging.
    field a("hcurl", spantree), v("h1"), x("x"), y("y"), z("z");
    
    // Use interpolation order 1:
    a.setorder(wholedomain, 1);
    v.setorder(wholedomain, 1);
    
    // Gauge the vector potential field on the whole volume:
    a.setgauge(wholedomain);
    
    // Put a magnetic wall (i.e. set field a to 0) all around the domain (no magnetic flux can cross it):
    a.setconstraint(domainskin);
    // Ground v on a node to have a well-defined problem:
    v.setconstraint(ground);
    
    
    // The magnetic vector potential a is rewritten as the sum of a known
    // source term and an unknown part to be determined: atotal = a + asource.
    // We want to apply an external induction field b of 1 mT in 
    // the z direction that linearly increases over time:
    expression bsource = array3x1(0,0,1e-3)*t();
    // Since b = curl(a) a valid choice for a is:
    expression asource = 1e-3*0.5*array3x1(-y,x,0)*t();
    // And dt(a)/dt is:
    expression dtasource = 1e-3*0.5*array3x1(-y,x,0);
    
    
    // The strong form of the magnetodynamic a-v formulation is 
    // 
    // curl( 1/mu * curl(a) ) + sigma * (dt(a) + grad(v)) = js, with b = curl(a) and e = -dt(a) - grad(v)
    //
    expression e = -dt(a)-dtasource - grad(v);
    
    // The conductivity of the high temperature superconductor is modeled using 
    // a power law relation between the electric field and the current density:
    //
    // J = Jc/Ec^(1/n) * norm(E)^((1-n)/n) * E
    //    |_______________________________|
    //                sigma(E)
    //
    // where the critical current density Jc [A/m^2] and critical electric field Ec [V/m] are supposed 
    // constant here and the J(E) relation is isotropic. An extra 1e-11 is added to never divide by 0.
    //
    double n = 30.0, Ec = 1e-4, Jc = 1e8;
    
    expression sigma = Jc/( pow(Ec,1.0/n) ) * pow( norm(e) + 1e-11, (1.0-n)/n );
    
    
    // Define the weak magnetodynamic formulation:
    formulation magdyn;
    
    // Magnetic equation (add an extra odd integration degree for convergence):
    magdyn += integral(wholedomain, 1/mu*( curl(dof(a)) + bsource ) * curl(tf(a)) );
    magdyn += integral(tube, sigma*( dt(dof(a)) + dtasource )*tf(a), +1 );
    magdyn += integral(tube, sigma*grad(dof(v))*tf(a) );
    // Electric equation:
    magdyn += integral(tube, sigma*grad(dof(v))*grad(tf(v)) );
    magdyn += integral(tube, sigma*( dt(dof(a)) + dtasource )*grad(tf(v)), +1 );
    
    
    // Start the implicit Euler time resolution with the field values and a zero time derivative:
    impliciteuler eul(magdyn, vec(magdyn));
    
    // Tolerance on the nonlinear iteration:
    eul.settolerance(1e-3);
    // Relaxation factor for the nonlinear iteration (decrease to improve convergence):
    eul.setrelaxationfactor(0.5);
    
    // Use a 5 seconds timestep and run until 150 sec (30 timesteps):
    int numsteps = 30;
    double timestep = 5.0;
    std::vector<double> bcenter(numsteps);
    
    settime(0);
    for (int i = 0; i < numsteps; i++)
    {
        // Use a fixed-point nonlinear iteration at every timestep (25 nl iterations max).
        eul.next(timestep, 25);
        
        // Write b = curl(a) to disk:
        norm(curl(a) + bsource).write(wholedomain, "normbtotal" + std::to_string(i+1000) + ".vtk"); 
        
        // Output the b induction field [T] at the tube center to assess the shielding effectiveness.
        bcenter[i] = norm(curl(a) + bsource).interpolate(wholedomain, {1e-10,0,0})[0];
        std::cout << std::endl << "b source " << timestep*(i+1) << " mT --> b tube center " << bcenter[i]*1e3 << " mT" << std::endl;
    }
    // Write to file the field at the tube center for all timesteps:
    writevector("bcenter.csv", bcenter);
    
    clk.print("Total computation time:");
    
    // Code validation line. Can be removed.
    std::cout << 1;
}

// THE MESH BELOW IS FULLY STRUCTURED AND IS CREATED USING THE (BASIC) SPARSELIZARD GEOMETRY CREATION TOOL.
// THE ADVANTAGE OF IT IS THAT THE CODE ABOVE CAN BE CALLED FOR ANY TUBE DIMENSION WITHOUT NEEDING CALLS TO EXTERNAL MESHING SOFTWARE.
// GMSH COULD HAVE BEEN USED AS AN ALTERNATIVE TO DEFINE THE GEOMETRY AND THE MESH.

mesh createmesh(double thtube, double htube, double rintube, double linf, int nhtube, int ncircumtube, int nthtube, int nair, int prog)
{
    int wholedomain = 1, tube = 2, insidetube = 3, tubeskin = 4, domainskin = 5, ground = 6;
    
    // To be sure to have nodes at z = 0:
    if (nhtube%2 == 0)
        nhtube++;
    
    field x("x"), y("y");
    
    // Create the footprint face for the air inside the tube:
    shape disk("disk", -1, {0,0,0}, rintube, ncircumtube*4);
    // Create the footprint face for the tube:
    shape arcin = disk.getsons()[0];
    shape arcout("arc", -1, {rintube+thtube,0,0, 0,rintube+thtube,0, 0,0,0}, ncircumtube+1);
    shape line1("line", -1, {rintube,0,0, rintube+thtube,0,0}, nthtube);
    shape line2 = line1.duplicate(); line2.rotate(0,0,90);
    
    std::vector<shape> quads(4);
    for (int i = 0; i < 4; i++)
    {
        shape curquad("quadrangle", -1, {arcin,line1,arcout,line2});
        curquad = curquad.duplicate();
        curquad.rotate(0,0,90.0*i);
        quads[i] = curquad;
    }
    shape tubefootprint("union", -1, quads);
    // Create the footprint for the remaining domain with a given progression 'prog':
    double rout = std::pow(linf*std::pow(rintube+thtube,prog-1), 1.0/prog);
    shape arcinf("arc", -1, {rout,0,0, 0,rout,0, 0,0,0}, ncircumtube+1);
    shape lineinf1("line", -1, {rintube+thtube,0,0, rout,0,0}, nair);
    shape lineinf2 = lineinf1.duplicate(); lineinf2.rotate(0,0,90);
    
    std::vector<shape> airquads(4);
    for (int i = 0; i < 4; i++)
    {
        shape curquad("quadrangle", -1, {arcout,lineinf1,arcinf,lineinf2});
        curquad = curquad.duplicate();
        curquad.rotate(0,0,90.0*i);
        airquads[i] = curquad;
    }
    shape airfootprint("union", -1, airquads);
    // Deform to have a mesh progression in the air:
    expression radius = sqrt(x*x+y*y);
    airfootprint.move(array3x1(x*pow(radius/(rintube+thtube),prog-1)-x,y*pow(radius/(rintube+thtube),prog-1)-y,0));
    
    
    // Create the cylinder region inside the tube:
    shape voltubeinside = disk.duplicate().extrude(insidetube, htube, nhtube);
    voltubeinside.shift(0,0,-htube/2);
    
    // Create the tube region:
    shape voltube = tubefootprint.duplicate().extrude(tube, htube, nhtube);
    voltube.shift(0,0,-htube/2);
    
    // Create the air cylinder around the tube:
    shape volairaroundtube = airfootprint.duplicate().extrude(-1, htube, nhtube);
    volairaroundtube.shift(0,0,-htube/2);
    
    // Create the air above and below the tube:
    shape quadall("union", -1, {disk,tubefootprint,airfootprint});
    shape volabove = quadall.duplicate().extrude(-1, linf, nair);
    volabove.shift(0,0,htube/2);
    shape volbelow = volabove.duplicate();
    volbelow.scale(1,1,-1);
    
    // Define the whole volume:
    shape wholevol("union", wholedomain, {voltube.duplicate(), voltubeinside.duplicate(), volairaroundtube.duplicate(), volabove.duplicate(), volbelow.duplicate()});
    
    mesh mymesh;
    
    mymesh.regionskin(tubeskin, tube);
    mymesh.regionskin(domainskin, wholedomain);
    mymesh.boxselection(ground, wholedomain, 0, {rintube+thtube-1e-10,rintube+thtube+1e-10, -1e-10,1e-10, -1e-10,1e-10});
    
    mymesh.load({voltube, voltubeinside, wholevol});
    
    return mymesh;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

