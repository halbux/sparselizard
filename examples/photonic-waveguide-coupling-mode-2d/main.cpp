// This code simulates the eigenmodes of two nearby rectangular SiN dielectric photonic 
// waveguides burried in a SiO2 clad. To do so the 2D cross sectional geometry is 
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
// The waveguides are 500 nm wide and 250 nm heigh. Their spacing is 100 nm.
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

#include "sparselizardbase.h"

using namespace mathop;

void sparselizard(void)
{	
    // Give names to the physical region numbers :
    int wavg1 = 1, wavg2 = 2, clad = 3, all = 4, bound = 5, wgskin = 6;

    mesh mymesh;
    mymesh.regionskin(bound, all);
    mymesh.load("optical_waveguide.msh");
    
    int waveguides = regionunion({wavg1, wavg2});
    
    // Waveguide boundary plotting:
    expression(1).write(regionintersection({waveguides, clad}), "skin.vtu");


    wallclock clk;

    // Edge shape functions 'hcurl' for the tranverse electric field Et.
    // Nodal shape functions 'h1' for the longitudinal electric field El.
    // Fields x, y and z are the x, y and z coordinate fields.
    field Et("hcurl"), El ("h1"), x("x"), y("y"),z("z");

    // Use interpolation order 2 on the whole domain:
    Et.setorder(all, 2);
    El.setorder(all, 2);

    //unit vectors used for projection along the X Y Z directions
    expression ex(3,1,{1,0,0});
    expression ey(3,1,{0,1,0});
    expression ez(3,1,{0,0,1});
    expression ex2(2,1,{1,0});
    expression ey2(2,1,{0,1});

    // material properties definition
    parameter n,epsr,mur;
    n|clad=1.4;
    n|wavg1=2;
    n|wavg2=2;
    epsr|all=n*n;
    mur|all=1;
    // light property
    double lmb = 680e-9, c = 299792458, f0=c/lmb,k0 = 2*getpi()/lmb;
    std::cout<<endl<<"lmb= "<< lmb*1e9<<" nm   f0= "<<f0<<" Hz   k0= "<<k0<<" m-1    "<<endl<<endl;

    // perfect conductor BC : nxE=0
    El.setconstraint(bound);
    Et.setconstraint(bound);

    // Weak formulation for the eigenvalue problem 
    // for confined propagation of an EM wave in a waveguide:
    //
    // 1/mur*curl(Et)*curl(tfEt)-k0^2*epsr*Et*tEt = -bt^2/mur*( grad(El)*tfEt) + Et*tfEt )
    //
    // bt^2/mur*( grad(El)*grad(tfEl) + Et*grad(tfEl) ) = k0^2*bt^2*epsr* El*tfEl
    //
    // the curl and grad operator are particylarized in the transverse plane, meaning X Y
    // further, the bt constant, coming from the space derivation of exp(i*Bt*z) with be replace by a time derivative dt
    //


    //expressions needed 
    expression dtdtgradEl(3,1,{dtdt(grad(dof(El)))*ex2,dtdt(grad(dof(El)))*ey2,0});
    expression gradtfEl(3,1,{grad(tf(El))*ex2,grad(tf(El))*ey2,0});

    formulation mode;

    mode += integral(all, curl(dof(Et))*curl(tf(Et))- k0*k0*mur*epsr*(dof(Et))*tf(Et));
    mode += integral(all, dtdtgradEl*tf(Et)+ dtdt(dof(Et))*tf(Et));

    mode += integral(all, dtdtgradEl*gradtfEl+dtdt(dof(Et))*gradtfEl);
    mode += integral(all, -k0*k0*mur*epsr*dtdt(dof(El))*tf(El));	

    // Generate the algebraic system Ax = b:
    mode.generate();

    // Get the stiffness, damping and mass matrix:
    mat K = mode.K();
    mat M = mode.M();

    // Remove the rows and columns corresponding to the 0 constraints:
    K.removeconstraints();
    M.removeconstraints();

    // Create the object to solve the eigenvalue problem 
    eigenvalue eig(K,M); 

    // Compute the 5 eigenvalues closest to the target magnitude bt^2 (propagation constant^2)
    // basically we look for mode around neff_guess, guess of the effective refarctive index
    double neff_guess=1.6, bt=k0*neff_guess, bt2=pow(bt,2);
    std::cout<<"Look around neff = "<<neff_guess<<"    beta = "<<bt<<" rad/m"<<endl<<endl;
    // propagation mode eigenvalue is purely imaginary : i*bt => look around -bt^2
    eig.compute(5, -bt2);

    // Catalogue eigenvectors and eigenvalues:
    std::vector<vec> myrealeigenvectors = eig.geteigenvectorrealpart();
    std::vector<vec> myimageigenvectors = eig.geteigenvectorimaginarypart();
    std::vector<double> myrealeigenvalues=eig.geteigenvaluerealpart();
    std::vector<double> myimageigenvalues=eig.geteigenvalueimaginarypart();

    int md=1;
    double btc;
    std::vector<double> neffcs;


    std::cout << "Modes retrieved:" <<endl<<endl;

    expression Etotal = array3x1(Et*ex, Et*ey, El);

    // Loop on all eigenvalues found:
    for (int i = 0; i < myrealeigenvalues.size(); i++)
    {
        if( myrealeigenvalues[i] < 0 )
        {
            // Transfer the data from the ith eigenvector to field Et and El:
            Et.setdata(all, myrealeigenvectors[i]);
            El.setdata(all, myrealeigenvectors[i]);

            //compute propagation constant and effective refractive index
            btc=sqrt(abs(myrealeigenvalues[i]));
            double neffc=btc/k0;
            // We need to separate into clad and waveguide regions to keep the discontinuity:
            Etotal.write(clad, "Eclad_"+ std::to_string(md) +".vtu", 2);
            Etotal.write(regionunion({wavg1,wavg2}), "Ewav_"+ std::to_string(md) +".vtu", 2);
            // display retrieved mode information
            std::cout<<"beta["<<md<<"] = "<<btc<<" rad/m    neff = "<<neffc<<endl;

            neffcs.push_back(neffc);
            
            md++;
        }
    }

    // Write a .pvd file to visualize all the modes together:
    grouptimesteps("Eclad.pvd", "Eclad_", 1, neffcs);
    grouptimesteps("Ewav.pvd", "Ewav_", 1, neffcs);

    std::cout<<endl;
    clk.print("Total time elapsed:");
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

