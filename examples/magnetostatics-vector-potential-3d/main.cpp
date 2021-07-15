// This code simulates the magnetic field created by an imposed current density
// in an aluminium wire when a magnetic shield (a steel cylinder) is placed nearby.
//
// This example illustrates the resolution of a general magnetostatic problem with
// current density sources on a 3D geometry. It uses the magnetic vector potential
// based on edge shape functions 'hcurl' and combined with a gauge. The gauge is
// needed because the algebraic system to solve is singular. The effect of the gauge
// is to remove all degrees of freedom associated to higher order gradient type shape 
// functions from the algebraic system. Additionally all degrees of freedom associated
// to lowest order (order 0) edge shape functions are removed for all edges in a
// spanning tree. After having removed these two types of degrees of freedom the 
// matrix becomes non-singular and can be solved with a usual direct solver.
// 
// It is worth noting that the spanning tree growth HAS TO START on the regions
// where the magnetic vector potential field 'a' is constrained. This is to avoid
// issues of forcing the value of the magnetic flux through faces where we do not 
// want it to be constrained.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
	// Give names to the region numbers defined in the mesh file:
	int conductor = 1, shield = 2, air = 3, contour = 4;

	// Load the mesh:
	mesh mymesh("magmesh.msh");

	// Define the whole volume region for convenience:
	int wholedomain = selectunion({conductor,shield,air});

	// Define the magnetic permeability mu [H/m] in the air, conductor (aluminium) and magnetic shield (steel):
	double mu0 = 4*getpi()*1e-7;
	
	parameter mu;

	mu|air = mu0;
	mu|conductor = mu0;
	mu|shield = 1000*mu0;

	// Define a spanning tree to gauge the magnetic vector potential (otherwise the matrix is singular).
	// Start growing the tree from the regions with constrained potential vector (here the contour): 
	spanningtree spantree({contour});
	// Write it for illustration:
	spantree.write("spantree.pos");

	// Define the magnetic vector potential 'a' and provide the tree. Use edge shape functions 'hcurl'.
	field a("hcurl", spantree);

	// Gauge the vector potential field on the whole volume:
	a.setgauge(wholedomain);

	// Use higher interpolation orders where needed:
	a.setorder(shield, 2);
	a.setorder(conductor, 1);
	a.setorder(air, 1);

	// Put a magnetic wall (i.e. set field 'a' to 0) all around the domain (no magnetic flux can cross it):
	a.setconstraint(contour);

	// Define the magnetostatic formulation:
	formulation magnetostatics;

	// The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
	magnetostatics += integral(wholedomain, 1/mu* curl(dof(a)) * curl(tf(a)) );
	// A current density of 1A/m^2 flows in the z direction in the conductor region:
	magnetostatics += integral(conductor, -array3x1(0,0,1)*tf(a));

	// Generate the algebraic matrices A and vector b of the Ax = b problem:
	magnetostatics.generate();

	// Get the x solution vector:
	vec sol = solve(magnetostatics.A(), magnetostatics.b());

	// Update field 'a' with the solution that has just been computed:
	a.setdata(wholedomain, sol);
	
	// Write the magnetic induction field b = curl(a) [T]:
	curl(a).write(wholedomain, "b.pos", 1);
	
	// Code validation line. Can be removed:
	std::cout << (norm(a).max(wholedomain,4)[0] < 1.96437e-08 && norm(a).max(wholedomain,4)[0] > 1.96435e-08);
}

