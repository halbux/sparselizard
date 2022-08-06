// This code simulates the AC magnetic field and induced currents created
// by three conducting wires, each with a prescribed total current flow.
//
// The so-called h-phi (a.k.a. t-w) formulation is used. The 'main' section
// of this code can be (and has been) used as it is for 3D h-phi simulations. 
//
// Credits: J. Ruuskanen


#include "sparselizard.h"
#include "gmsh.h"


using namespace sl;

int air = 1, cond1 = 2, cond2 = 3, cond3 = 4, cohomcut1 = 5, cohomcut2 = 6, cohomcut3 = 7, point = 8;

mesh createmesh(void);

int main(void)
{
    // Create the mesh with the gmsh api:
    mesh mymesh = createmesh();

    printversion();

    int all = selectall();
    int conductor = selectunion({cond1, cond2, cond3});
    int cohomcuts = selectunion({cohomcut1, cohomcut2, cohomcut3});

    field x("x");
    
    double mu0 = 4*getpi()*1e-7;

    // Define the electric conductivity and magnetic permeability:
    parameter sigma, mu;

    // Coordinate dependent conductivity for illustration:
    sigma|conductor = 1e7*(abs(x+0.014)/0.01);

    mu|air = mu0;
    mu|conductor = 1000*mu0;
    
    sigma.write(conductor, "sigma.pos", 1);

    // Set the working frequency to 50 Hz:
    setfundamentalfrequency(50);

    int condskin = selectintersection({air, conductor}, 1);

    // Since the solution has a component in phase with the actuation
    // and a quadrature component we need 2 harmonics at 50Hz 
    // (harmonic 1 is DC, 2 is sine at 50Hz and 3 cosine at 50Hz).
    std::vector<int> harms = {2,3};

    field hc("hcurl", harms);

    hc.setorder(conductor, 2);
    hc.setconstraint(condskin);

    // Source magnetic field:
    field hs("hcurl", harms);

    hs.setorder(all, 0);

    // Set the sources [A] on the three cohomology cuts for the 50Hz sine harmonic.
    //
    // The total current in each conductor and the cohomology sources are generally not equal.
    // The actual total current flowing in each conductor is printed at the end of this code.
    // It is determined by the circulation on a closed loop around each conductor. Since the
    // orientation of a cohomology cut depends on the mesh, it is recommended to visualize
    // the source field 'hs.pos' to determine which cohomology source coefficients must be
    // used to get the desired total current in each conductor. 
    //
    hs.harmonic(2).setcohomologysources({cohomcut1, cohomcut2, cohomcut3}, {3, 1, -2});

    hs.write(cohomcuts, "hs.pos", 1);

    // Scalar potential of the h-phi formulation:
    field phi("h1", harms);

    phi.setorder(all, 2);   
    // Fix the potential to 0 at a point in space:
    phi.setconstraint(point);

    // The gradient has only two components in 2D but 'hcurl' fields have three.
    // Add a 0-valued component to have a three component gradient.
    // 
    // In the formulation below the contribution of phi to the conductor region elements
    // should be taken into account but there should be no phi unknown associated to the
    // conductor except for its skin region. To achieve this an extra argument is
    // provided to the dof() and tf() functions to define them only on the air region.
    //
    expression graddofphi = grad(dof(phi, air)).resize(3,1);
    expression gradtfphi = grad(tf(phi, air)).resize(3,1);

    formulation magdyn;

    // The strong form of the magnetodynamic h-phi formulation is
    //
    // curl( 1/sigma * curl(h) ) + dt(b) = 0
    //
    // with hair = hs + grad(phi) and hconductor = hc + hs + grad(phi)
    //
    magdyn += integral(conductor, 1/sigma * (curl(dof(hc)) + curl(hs)) * curl(tf(hc)));

    magdyn += integral(conductor, mu * (dt(dof(hc)) + dt(hs) + dt(graddofphi)) * (tf(hc) + gradtfphi));
    magdyn += integral(air, mu * (dt(hs) + dt(graddofphi)) * gradtfphi);
    // Avoid 0 diagonal entries in the assembled matrix:
    magdyn += integral(air, 1e-50 * graddofphi * gradtfphi);

    magdyn.solve();

    expression gradphi = grad(phi).resize(3,1);

    // Current density:
    expression j = curl(hc) + curl(hs);
    // Air and conductor magnetic field:
    expression hair = hs + gradphi;
    expression hcond = hc + hs + gradphi;

    j.write(conductor, "j.pos", 2);

    (mu*hcond).write(conductor, "bcond.pos", 2);
    (mu*hair).write(air, "bair.pos", 2);
    
    double I1 = getharmonic(2, compz(j)).integrate(cond1, 5);
    double I2 = getharmonic(2, compz(j)).integrate(cond2, 5);
    double I3 = getharmonic(2, compz(j)).integrate(cond3, 5);
    
    std::cout << "Total current in wire 1/2/3 is " << I1 << "/" << I2 << "/" << I3 << " A" << std::endl;

    // It can be seen below that the h-phi formulation gives a better accuracy on b than on j.
    double B2max = norm(getharmonic(2, mu*hcond)).max(conductor, 5)[0];
    double B3max = norm(getharmonic(3, mu*hcond)).max(conductor, 5)[0];
    std::cout << "B in-phase/quadrature max is " << B2max << " / " << B3max << " T" << std::endl;
    
    double jz2max = norm(getharmonic(2, compz(j))).max(conductor, 5)[0];
    double jz3max = norm(getharmonic(3, compz(j))).max(conductor, 5)[0];
    std::cout << "J in-phase/quadrature max is " << jz2max << " / " << jz3max << " A/m2" << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (std::abs(I1-3)/3 < 1e-12 && std::abs(I2+4)/4 < 1e-12 && std::abs(I3-2)/2 < 1e-12 && std::abs(B2max-0.4303)/0.4303 < 1e-3 && std::abs(B3max-0.1533)/0.1533 < 1e-3 && std::abs(jz2max-690000)/690000 < 1e-3 && std::abs(jz3max-576000)/576000 < 1e-3);
}

mesh createmesh(void)
{
    // Mesh size:
    double refinefact = 8.0;
    double msair = 1e-2 / refinefact, mscond = 5e-4 / refinefact;

    // Radii of the air domain and the conductors:
    double rair = 5e-2, rcond = 2e-3;
    
    // Conductor positions:
    std::vector<std::vector<double>> xyzc = {{-rcond/2-5*rcond, -rcond/2, 0}, {-rcond/2, -rcond/2, 0}, {-rcond/2+5*rcond, -rcond/2, 0}};

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 2);

    std::vector<std::pair<int, int>> ov;
    std::vector<std::vector<std::pair<int, int>>> ovv;

    // Define the geometry:
    int airgeo = gmsh::model::occ::addDisk(0, 0, 0, rair, rair);
    int cond1geo = gmsh::model::occ::addDisk(xyzc[0][0], xyzc[0][1], xyzc[0][2], rcond, rcond);
    int cond2geo = gmsh::model::occ::addDisk(xyzc[1][0], xyzc[1][1], xyzc[1][2], rcond, rcond);
    int cond3geo = gmsh::model::occ::addDisk(xyzc[2][0], xyzc[2][1], xyzc[2][2], rcond, rcond);

    gmsh::model::occ::fragment({{2, airgeo}}, {{2, cond1geo}, {2, cond2geo}, {2, cond3geo}}, ov, ovv);

    gmsh::model::occ::synchronize();

    airgeo = get<1>(ov[3]);
    cond1geo = get<1>(ov[0]);
    cond2geo = get<1>(ov[1]);
    cond3geo = get<1>(ov[2]);

    gmsh::model::occ::synchronize();

    gmsh::model::addPhysicalGroup(2, {airgeo}, air);
    gmsh::model::addPhysicalGroup(2, {cond1geo}, cond1);
    gmsh::model::addPhysicalGroup(2, {cond2geo}, cond2);
    gmsh::model::addPhysicalGroup(2, {cond3geo}, cond3);

    // Set mesh size:
    std::vector<std::pair<int, int> > domainbnd, condbnd;
    
    gmsh::model::getBoundary({{2, airgeo}}, domainbnd, true, true, true);
    gmsh::model::mesh::setSize(domainbnd, msair);

    gmsh::model::getBoundary({ {2,cond1geo} }, condbnd, true, true, true);
    gmsh::model::mesh::setSize(condbnd, mscond);

    gmsh::model::getBoundary({ {2,cond2geo} }, condbnd, true, true, true);
    gmsh::model::mesh::setSize(condbnd, mscond);

    gmsh::model::getBoundary({ {2,cond3geo} }, condbnd, true, true, true);
    gmsh::model::mesh::setSize(condbnd, mscond);

    // Compute the cohomology basis:
    gmsh::model::mesh::addHomologyRequest("Cohomology", {air}, {}, {1});

    // Mesh in 2D:
    gmsh::model::mesh::generate(2);

    // Load mesh in sparselizard:
    mesh mymesh;
    mymesh.selectanynode(point, air);
    mymesh.load("gmsh:api", 1);

    gmsh::finalize();

    mymesh.write("wires2d.msh");

    std::cout << "Mesh refinement factor is " << refinefact << "." << std::endl;

    return mymesh;
}

