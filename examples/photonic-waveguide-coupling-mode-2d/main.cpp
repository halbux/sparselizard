// This code simulates the eigenmodes of two nearby rectangular SiN dielectric photonic 
// waveguides buried in a SiO2 clad. To do so the 2D cross sectional geometry is 
// considered and the propagation is assumed to be along the z axis.
// To solve this eigenvalue problem in 2D while taking into account the z component
// the electric field is split into a transverse and a longitudinal part.
// Since the searched modes are confined in and around the waveguide, a perfect 
// conductor boundary condition is valid, as long as the domain is large compared 
// to the waveguide dimensions.
//
// More details can be found in 'NASA technical paper 3485'.
//
//
// The waveguides are 500 nm wide and 250 nm high. Their spacing is 100 nm.
//       _____________________________________________________
//      |                                                     |
//      |   SiO2      ___________     ___________             |
//      |            |           |   |           |            |
//      |            |    SiN    |   |    SiN    |            |
//      |            |___________|   |___________|            |
//      |                                                     |
//      |_____________________________________________________|
//
// Credits: R. Haouari

#include "sparselizard.h"

using namespace sl;

int main(void)
{	
    // Give names to the physical region numbers :
    int wavg1 = 1, wavg2 = 2, clad = 3, all = 4, bound = 5;

    mesh mymesh;
    mymesh.selectskin(bound, all);
    mymesh.load("optical_waveguide.msh");
    
    int waveguides = selectunion({wavg1, wavg2});
    
    // Waveguide boundary plotting:
    expression(1).write(selectintersection({waveguides, clad}, 1), "skin.vtu", 1);

    wallclock clk;

    // Edge shape functions 'hcurl' for the tranverse electric field Et.
    // Nodal shape functions 'h1' for the longitudinal electric field El.
    field Et("hcurl"), El("h1");

    // Use interpolation order 2 on the whole domain:
    Et.setorder(all, 2);
    El.setorder(all, 2);

    // Material properties definition
    parameter n, epsr, mur;
    n|clad = 1.4;
    n|wavg1 = 2.0;
    n|wavg2 = 2.0;
    epsr|all = n*n;
    mur|all = 1.0;
    
    // Light property
    double lambda = 680e-9, c = 299792458, f0 = c/lambda, k0 = 2.0*getpi()/lambda;
    std::cout << std::endl << "lambda = " << lambda*1e9 <<" nm, f0 = " << f0 << " Hz, k0 = " << k0 << " 1/m" << std::endl << std::endl;

    // Perfect conductor boundary condition:
    El.setconstraint(bound);
    Et.setconstraint(bound);

    // Weak formulation for the eigenvalue problem for confined propagation of an EM wave in a waveguide:
    //
    // 1/mur*curl(Et)*curl(Et') - k0^2*epsr*Et*Et' = -bt^2/mur*( grad(El)*Et') + Et*Et' )
    //
    // bt^2/mur*( grad(El)*grad(El') + Et*grad(El') ) = k0^2*bt^2*epsr* El*El'
    //
    // where curl() and grad() have a special definition in the transverse (xy) plane.
    //
    // The bt constant coming from the space derivation of exp(i*bt*z) is replaced by a time derivative dt().


    // Operators grad() and curl() in the transverse plane:
    expression dtdtgradEl(3,1,{dtdt(dx(dof(El))), dtdt(dy(dof(El))), 0});
    expression gradtfEl(3,1,{dx(tf(El)), dy(tf(El)), 0});

    formulation mode;

    mode += integral(all, curl(dof(Et))*curl(tf(Et))- k0*k0*mur*epsr*(dof(Et))*tf(Et));
    mode += integral(all, dtdtgradEl*tf(Et)+ dtdt(dof(Et))*tf(Et));

    mode += integral(all, dtdtgradEl*gradtfEl+dtdt(dof(Et))*gradtfEl);
    mode += integral(all, -k0*k0*mur*epsr*dtdt(dof(El))*tf(El));	

    mode.generate();

    // Get the stiffness matrix K and mass matrix M:
    mat K = mode.K();
    mat M = mode.M();

    // Create the object to solve the eigenvalue problem:
    eigenvalue eig(K, M); 

    // Compute the 5 eigenvalues closest to the target magnitude bt^2 (propagation constant^2).
    // We are looking for modes around neff_target, the target effective refractive index.
    double neff_target = 1.6, bt = k0*neff_target, bt2 = std::pow(bt,2.0);
    std::cout << "Target is neff = " << neff_target << ", beta = " << bt << " rad/m" << std::endl << std::endl;
    // Propagation mode eigenvalue is purely imaginary (i*bt), we thus target -bt^2:
    eig.compute(5, -bt2);

    // Get all eigenvectors and eigenvalues found:
    std::vector<vec> myrealeigenvectors = eig.geteigenvectorrealpart();
    std::vector<double> myrealeigenvalues = eig.geteigenvaluerealpart();


    std::vector<double> neffcs;
    
    expression Etotal = array3x1(compx(Et), compy(Et), El);

    // Loop on all eigenvalues found:
    int index = 1;
    for (int i = 0; i < myrealeigenvalues.size(); i++)
    {
        if (myrealeigenvalues[i] < 0)
        {
            // Transfer the data from the ith eigenvector to fields Et and El:
            Et.setdata(all, myrealeigenvectors[i]);
            El.setdata(all, myrealeigenvectors[i]);

            // Compute the propagation constant and the effective refractive index:
            double btc = std::sqrt(std::abs(myrealeigenvalues[i]));
            double neffc = btc/k0;
            // We need to write separately on the clad and waveguide regions to visualize the discontinuity:
            Etotal.write(clad, "Eclad_"+ std::to_string(index) +".vtu", 2);
            Etotal.write(selectunion({wavg1,wavg2}), "Ewav_"+ std::to_string(index) +".vtu", 2);
            // Display mode information:
            std::cout << "Mode " << index << ": beta = " << btc << " rad/m, neff = " << neffc << std::endl;

            neffcs.push_back(neffc);
            
            index++;
        }
    }

    // Write a .pvd file to visualize all the modes together:
    grouptimesteps("Eclad.pvd", "Eclad_", 1, neffcs);
    grouptimesteps("Ewav.pvd", "Ewav_", 1, neffcs);

    std::cout << std::endl;
    clk.print("Total time elapsed:");

    // Code validation line. Can be removed.
    std::cout << (neffcs[0] < 1.65780 && neffcs[0] > 1.65778);
}

