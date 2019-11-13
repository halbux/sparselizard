// This code simulates the coupling between 2 dielectric photonic waveguide.
// To do so, we need to evaluate the propagation modes of the waveguides whoose 
// cross sectional geometry is represented in 2D.
// Propagation is supposed to be along the z axis.
// The geometry considered is a rectangular waveguide burried in a clad.
// An optically isotropic material is considered here.
// To solve this eigenvalue problem and lift undtermination regarding the Z component,
// a split of the electric field into transverses and longitudinal components is needed.
// Since the desired mode is supposed to be confined in the wageguide
// a perfect conductor BC type set at the outer boundary is valid, as long as 
// the domain is big enough compared to the waveguide dimensions.
//
// Credits: R. Haouari

#include "sparselizardbase.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace mathop;

// A gives the waveguide width
// B gives the waveguide height
// D is the distance bewteen waveguides
mesh createmesh(double a, double b, double d);

void sparselizard(void)
{	
    // Give names to the physical region numbers :
    int wavg1 = 1, wavg2=2, cross = 3,bound=4, wavg_bound=5;
    //Create mesh
    mesh mymesh = createmesh(.5e-6,.25e-6,100e-9);
    mymesh.write("opticalwaveguide.msh");
    //Create new regions through boolean operation:
    int clad=regionexclusion(cross, regionunion({wavg1,wavg2})), edges=regionunion({bound,wavg_bound});
    //recover the edges lost from the regionexclusion process 
    clad=regionunion({clad,wavg_bound});
    //geometry edge plotting
    expression(1).write(edges,"edge.vtu");


    wallclock clk;

    // Edge shape functions 'hcurl' for the tranverse electric field Et.
    // Nodal shape functions 'h1' for the longitudinal electric field El.
    // Fields x, y and z are the x, y and z coordinate fields.
    field Et("hcurl"), El ("h1"), x("x"), y("y"),z("z"),nor("h1xyz");

    // Use interpolation order 2 on the whole domain:
    Et.setorder(cross, 2);
    El.setorder(cross, 2);
    nor.setorder(bound,2);

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
    epsr|cross=n*n;
    mur|cross=1;
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

    mode += integral(cross, curl(dof(Et))*curl(tf(Et))- k0*k0*mur*epsr*(dof(Et))*tf(Et));
    mode += integral(cross, dtdtgradEl*tf(Et)+ dtdt(dof(Et))*tf(Et));

    mode += integral(cross, dtdtgradEl*gradtfEl+dtdt(dof(Et))*gradtfEl);
    mode += integral(cross, -k0*k0*mur*epsr*dtdt(dof(El))*tf(El));	

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
            Et.setdata(cross, myrealeigenvectors[i]);
            El.setdata(cross, myrealeigenvectors[i]);

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

mesh createmesh(double a, double b, double d)
{
    // Give names to the physical region numbers:
    int wavg1 = 1, wavg2=2, cross = 3,bound=4, wavg_bound=5;

    int sclx=7,scly=7;//choose odd integer number !!!
    double dx=sclx*a, dy=scly*b;

    // Number of mesh layers in the waveguide:
    int nx = 21,ny=21;


    // whole cross-section
    shape whole("quadrangle", cross , {-dx/2,-dy/2,0, dx/2,-dy/2,0, dx/2,dy/2,0, -dx/2,dy/2,0}, {sclx*(nx-1)+1,scly*(ny-1)+1,sclx*(nx-1)+1,scly*(ny-1)+1});
    shape guide1("quadrangle", wavg1 , {-a-d/2,-b/2,0, -d/2,-b/2,0, -d/2,b/2,0, -a-d/2,b/2,0}, {nx,ny,nx,ny});
    shape guide2("quadrangle", wavg2 , {d/2,-b/2,0, a+d/2,-b/2,0, a+d/2,b/2,0, d/2,b/2,0}, {nx,ny,nx,ny});


    // Provide to the mesh all shapes of interest:
    mesh mymesh;

    mymesh.regionskin(bound,cross);
    mymesh.regionskin(wavg_bound,wavg1);
    mymesh.regionskin(wavg_bound,wavg2);

    mymesh.load({whole,guide1,guide2});

    return mymesh;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    sparselizard();

    SlepcFinalize();

    return 0;
}

