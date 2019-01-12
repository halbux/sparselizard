// This code simulates the static magnetic field created by 5 permanent magnets (1cm x 1cm) placed side by side.
// A steel region above the magnets is perturbing the magnetic field lines.
//
// The magnetic scalar potential is used to solve this problem. This is valid since there are no current sources.
// The permanent magnets are treated as pre-magnetised pieces of non-magnetic material (since all magnetic domains are already fully oriented).


#include "sparselizardbase.h"


using namespace mathop;

void sparselizard(void)
{	
	// The domain regions as defined in 'halbacharray.geo':
	int magnet1 = 1, magnet2 = 2, magnet3 = 3, magnet4 = 4, magnet5 = 5, steel = 6, air = 7, zeropotential = 8;

	// The mesh can be curved!
	mesh mymesh("halbacharray.msh");

	// Define new physical regions for convenience:
	int magnets = regionunion({magnet1, magnet2, magnet3, magnet4, magnet5});
	int wholedomain = regionunion({magnets, steel, air});

	// Nodal shape functions 'h1' for the magnetic scalar potential 'phi'.
	field phi("h1");

	// Use interpolation order 2 on the whole domain:
	phi.setorder(wholedomain, 2);

	// The magnetic scalar potential (just like the electric potential)
	// needs to be fixed at least at one node to have a well-posed problem.
	// Here it is forced to 0 at a permanent magnet corner point.
	phi.setconstraint(zeropotential);

	// Vacuum magnetic permeability [H/m]: 
	double mu0 = 4*getpi()*1e-7;

	// Define the permeability in all regions:
	parameter mu;

	mu|air = mu0;
	mu|steel = 1000*mu0;
	mu|magnets = mu0;


	formulation magnetostatic;

	// In the absence of current sources the static Maxwell equations give:
	// 
	// curl h = 0
	// 
	// One can thus define a magnetic scalar potential 'phi' such that h = -grad(phi)
	//
	// Considering also that div(b) = 0 we get
	//
	// div(mur*(-grad(phi))) = 0
	//
	// with b = mur * h.
	//
	// In the permanent magnet region b = mu0 * (h + m),
	// i.e. the material is non-magnetic but it is pre-magnetised by the magnetisation vector m [A/m].
	// We thus get:
	// div(mu0*(-grad(phi)) + mu0*m) = 0
	// 
	// For a permanent magnet with a 1 Tesla magnetic induction b in the y direction we have:
	//
	// m = [0, 1 Tesla / mu0] or about [0, 800e3] A/m

	// The weak form corresponding to the above equations:
	magnetostatic += integral(wholedomain, -grad(dof(phi)) * mu * grad(tf(phi)) );

	// This is when all magnets are oriented in the y direction:
	magnetostatic += integral(magnets, array2x1(0, 800e3) * mu * grad(tf(phi)) );
	// This is in Halbach configuration (to maximise the magnetic field above the array):
	//magnetostatic += integral(magnet1, array2x1(-800e3, 0) * mu * grad(tf(phi)) );
	//magnetostatic += integral(magnet2, array2x1(0, -800e3) * mu * grad(tf(phi)) );
	//magnetostatic += integral(magnet3, array2x1(800e3, 0) * mu * grad(tf(phi)) );
	//magnetostatic += integral(magnet4, array2x1(0, 800e3) * mu * grad(tf(phi)) );
	//magnetostatic += integral(magnet5, array2x1(-800e3, 0) * mu * grad(tf(phi)) );


	magnetostatic.generate();

	// Solve the algebraic system Ax = b:
	vec sol = solve(magnetostatic.A(), magnetostatic.b());

	// Transfer the data from the solution vector to the phi field:
	phi.setdata(wholedomain, sol);

	// Write the magnetic scalar potential and the magnetic field with an order 2 interpolation.
	phi.write(wholedomain, "phi.pos", 2);
	norm(-grad(phi)).write(wholedomain, "hnorm.pos", 2);
	(-grad(phi)).write(wholedomain, "h.pos");
		
	// Evaluate the magnetic field 1cm above the center of the magnet array:
	std::vector<double> magfieldnorm = norm(-grad(phi)).interpolate(wholedomain, {0,0.02,0});
	std::cout << "Magnetic field 1cm above the array center: " << magfieldnorm[0] << " A/m" << std::endl << std::endl;

	// Write 20 magnetic field lines starting on the top side of the magnet array:
	shape ln("line", -1, {-0.025,0.01,0, 0.025,0.01,0}, 20);
	(-grad(phi)).streamline(regionexclusion(wholedomain, magnets), "magneticfieldline.pos", ln.getcoords(), 0.01/5);

	// Code validation line. Can be removed.
	std::cout << (magfieldnorm[0] < 64963.8 && magfieldnorm[0] > 64963.6);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

