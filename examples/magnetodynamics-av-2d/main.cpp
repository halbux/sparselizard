// This code simulates the AC magnetic field and induced currents created
// by three conducting wires, each with a prescribed total current flow.
//
// The so-called a-v (a.k.a. a-phi) formulation is used.


#include "sparselizard.h"
#include "gmsh.h"


using namespace sl;

int air = 1, cond1 = 2, cond2 = 3, cond3 = 4, skin = 9;

mesh createmesh(void);

int main(void)
{
    // Create the mesh with the gmsh api:
    mesh mymesh = createmesh();

    printversion();

    int all = selectall();
    int conductor = selectunion({cond1, cond2, cond3});

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

    // Since the solution has a component in phase with the actuation
    // and a quadrature component we need 2 harmonics at 50Hz 
    // (harmonic 1 is DC, 2 is sine at 50Hz and 3 cosine at 50Hz).
    std::vector<int> harms = {2,3};
    
    // Nodal shape functions 'h1' for the z component of the vector potential:
    field az("h1", harms);
    
    // In 2D the vector potential only has a z component:
    expression a = array3x1(0,0,az);

    // Set the interpolation order for field az:
    az.setorder(all, 2);
    
    // Put a magnetic wall (i.e. set field az to 0) all around the domain (no magnetic flux can cross it):
    az.setconstraint(skin);
    
    // In 2D the electric potential field v is not needed but in order to set a total current
    // flow condition on each conductor a field 'gradvz' equal to compz(grad(v)) is defined.
    field gradvz("h1", harms);
    
    // A port will be defined on each conductor so the lowest order is enough:
    gradvz.setorder(conductor, 1);
    
    // One primal/dual port pair per conductor. The primal is the lumped gradvz field on a region while the dual is
    // the total contribution of the Neumann term on that region, here the total current [A] in the z direction.
    port gradvz1(harms), Iz1(harms), gradvz2(harms), Iz2(harms), gradvz3(harms), Iz3(harms);
    
    gradvz.setport(cond1, gradvz1, Iz1);
    gradvz.setport(cond2, gradvz2, Iz2);
    gradvz.setport(cond3, gradvz3, Iz3);
    
    // Define the weak magnetodynamic formulation with a total 3 / -4 / 2 [A] flowing through conductor 1 / 2 / 3. 
    formulation magdyn;
    
    magdyn += Iz1.harmonic(2) - 3.0;
    magdyn += Iz1.harmonic(3) - 0;
    magdyn += Iz2.harmonic(2) + 4.0;
    magdyn += Iz2.harmonic(3) - 0;
    magdyn += Iz3.harmonic(2) - 2.0;
    magdyn += Iz3.harmonic(3) - 0;
    
    // The strong form of the magnetodynamic a-v formulation is 
    // 
    // curl( 1/mu * curl(a) ) + sigma * (dt(a) + grad(v)) = js
    //
    // with b = curl(a) and e = -dt(a) - grad(v)
    //
    // Magnetic equation:
    magdyn += integral(all, 1/mu * curl(dof(a)) * curl(tf(a)) );
    magdyn += integral(conductor, sigma * dt(dof(a)) * tf(a) + sigma * dof(gradvz) * tf(az) );
    // Electric equation:
    magdyn += integral(conductor, sigma * dof(gradvz) * tf(gradvz) + sigma * dt(dof(az)) * tf(gradvz) );
    
    magdyn.solve();
    
    expression ez = -dt(az) - gradvz;
    expression jz = sigma*ez;
    
    curl(a).write(conductor, "bcond.pos", 2);
    curl(a).write(air, "bair.pos", 2);
    array3x1(0,0,jz).write(conductor, "j.pos", 2);

    // Confirm the total current flow in each conductor:
    double I1 = getharmonic(2, jz).integrate(cond1, 5);
    double I2 = getharmonic(2, jz).integrate(cond2, 5);
    double I3 = getharmonic(2, jz).integrate(cond3, 5);
    
    std::cout << "Total current in wire 1/2/3 is " << I1 << "/" << I2 << "/" << I3 << " A" << std::endl;

    // It can be seen below that the a-v formulation gives a better accuracy on j than on b.
    double B2max = norm(getharmonic(2, curl(a))).max(conductor, 5)[0];
    double B3max = norm(getharmonic(3, curl(a))).max(conductor, 5)[0];
    std::cout << "B in-phase/quadrature max is " << B2max << " / " << B3max << " T" << std::endl;

    double jz2max = norm(getharmonic(2, jz)).max(conductor, 5)[0];
    double jz3max = norm(getharmonic(3, jz)).max(conductor, 5)[0];
    std::cout << "J in-phase/quadrature max is " << jz2max << " / " << jz3max << " A/m2" << std::endl;
    
    // Code validation line. Can be removed.
    std::cout << (std::abs(I1-3)/3 < 1e-12 && std::abs(I2+4)/4 < 1e-12 && std::abs(I3-2)/2 < 1e-12 && std::abs(B2max-0.430244)/0.430244 < 6e-3 && std::abs(B3max-0.153331)/0.153331 < 2e-3 && std::abs(jz2max-688904)/688904 < 9e-5 && std::abs(jz3max-576794)/576794 < 8e-5);
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
    std::vector<std::pair<int, int> > domainbnd, ironbnd;
    
    gmsh::model::getBoundary({{2, airgeo}}, domainbnd, true, true, true);
    gmsh::model::mesh::setSize(domainbnd, msair);

    gmsh::model::getBoundary({ {2,cond1geo} }, ironbnd, true, true, true);
    gmsh::model::mesh::setSize(ironbnd, mscond);

    gmsh::model::getBoundary({ {2,cond2geo} }, ironbnd, true, true, true);
    gmsh::model::mesh::setSize(ironbnd, mscond);

    gmsh::model::getBoundary({ {2,cond3geo} }, ironbnd, true, true, true);
    gmsh::model::mesh::setSize(ironbnd, mscond);

    // Mesh in 2D:
    gmsh::model::mesh::generate(2);

    // Load mesh in sparselizard:
    mesh mymesh;
    mymesh.selectskin(skin);
    mymesh.load("gmsh:api", 1);

    gmsh::finalize();

    mymesh.write("wires2d.msh");
    
    std::cout << "Mesh refinement factor is " << refinefact << "." << std::endl;

    return mymesh;
}

