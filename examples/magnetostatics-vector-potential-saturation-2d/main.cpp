// This code illustrates how to use a Newton iteration in the case of an anhysteretic magnetic
// saturation simulation using the a-v formulation on a simple 2D example. The 3D case can be
// obtained with minor adaptations to this code.
// An explicit implementation of the Newton iteration is proposed by default. An alternative
// implementation working with the fields instead of their increment is also mentionned. The
// latter leads to a simple nonlinear resolution with a loop on a regular .solve() call. 


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions as defined in 'twoconductors2d.geo':
    int conductor1 = 1, conductor2 = 2, insulator = 3, steel = 4, air = 5, bnd = 6;
    
    // Load the mesh while defining the domain boundary line:
    mesh mymesh;
    mymesh.selectskin(bnd);
    mymesh.load("magsat2d.msh");
    
    int wholedomain = selectall();
    int notmagnetic = selectunion({air, conductor1, conductor2,insulator});
    
    // The magnetic vector potential needs only a z component in 2D:
    field az("h1");
    az.setorder(wholedomain, 2);
    
    // Set a magnetic wall at the domain boundary:
    az.setconstraint(bnd);
    
    // Magnetic vector potential 'a' and magnetic induction 'b':
    expression a = array3x1(0,0,az);
    expression b = curl(a);
    
    // Measured data for the magnetic field h [A/m] and for b [T]:
    std::vector<double> hdata = {
    0.0000e+00, 5.5023e+00, 1.1018e+01, 1.6562e+01, 2.2149e+01, 2.7798e+01, 3.3528e+01,
    3.9363e+01, 4.5335e+01, 5.1479e+01, 5.7842e+01, 6.4481e+01, 7.1470e+01, 7.8906e+01,
    8.6910e+01, 9.5644e+01, 1.0532e+02, 1.1620e+02, 1.2868e+02, 1.4322e+02, 1.6050e+02,
    1.8139e+02, 2.0711e+02, 2.3932e+02, 2.8028e+02, 3.3314e+02, 4.0231e+02, 4.9395e+02,
    6.1678e+02, 7.8320e+02, 1.0110e+03, 1.3257e+03, 1.7645e+03, 2.3819e+03, 3.2578e+03,
    4.5110e+03, 6.3187e+03, 8.9478e+03, 1.2802e+04, 1.8500e+04, 2.6989e+04, 3.9739e+04,
    5.9047e+04, 8.8520e+04, 1.3388e+05, 2.0425e+05, 3.1434e+05, 4.8796e+05, 7.6403e+05};
    std::vector<double> bdata = {
    0.0000e+00, 5.0000e-02, 1.0000e-01, 1.5000e-01, 2.0000e-01, 2.5000e-01, 3.0000e-01,
    3.5000e-01, 4.0000e-01, 4.5000e-01, 5.0000e-01, 5.5000e-01, 6.0000e-01, 6.5000e-01,
    7.0000e-01, 7.5000e-01, 8.0000e-01, 8.5000e-01, 9.0000e-01, 9.5000e-01, 1.0000e+00,
    1.0500e+00, 1.1000e+00, 1.1500e+00, 1.2000e+00, 1.2500e+00, 1.3000e+00, 1.3500e+00,
    1.4000e+00, 1.4500e+00, 1.5000e+00, 1.5500e+00, 1.6000e+00, 1.6500e+00, 1.7000e+00,
    1.7500e+00, 1.8000e+00, 1.8500e+00, 1.9000e+00, 1.9500e+00, 2.0000e+00, 2.0500e+00,
    2.1000e+00, 2.1500e+00, 2.2000e+00, 2.2500e+00, 2.3000e+00, 2.3500e+00, 2.4000e+00};

    // Use a cubic spline interpolation for h(b):
    spline hbcurve(bdata, hdata);
    
    // Define nu as h = nu * b and fix the 0/0 division for the first entry:
    std::vector<double> nudata(hdata.size());
    for (int i = 0; i < hdata.size(); i++)
        nudata[i] = hdata[i]/bdata[i];
    nudata[0] = nudata[1];
    // Use a cubic spline interpolation for nu(b):
    spline nucurve(bdata, nudata);
    
    // Define nu on each region:
    parameter nu;

    nu|notmagnetic = 1/(4*getpi()*1e-7);
    nu|steel = expression(nucurve, norm(b)); // norm(b) is the abscissa
    
    expression h = nu * b;
    expression dhdb(hbcurve.getderivative(), norm(b));
    // Visualize how dhdb is interpolated using cubic splines:
    hbcurve.getderivative().write("dhdb.txt", 5);
    

    formulation magnetostatics;
    
    // Maxwell's equations in static form state that div(b) = 0 and curl(h) = j.
    // Equation div(b) = 0 allows to define a vector potential a such that b = curl(a).
    // Linearizing the nonlinear h(b) relation at b0 gives h = h(b0) + dh(b0)/db * db so
    // that equation curl(h) = j can be rewritten as curl(h(b0) + dh(b0)/db * db) = j,
    // which is the strong form for this simulation. In the non-magnetic regions we
    // have h = 0 + nu * (b - 0) = nu * b and thus the classical linear magnetostatic
    // strong form curl(nu * curl(a)) = j is obtained.
    
    // The following two terms can be used to avoid working with a Newton increment:
    // magnetostatics += integral(notmagnetic, nu * curl(dof(a)) * curl(tf(a)));
    // magnetostatics += integral(steel, (h + dhdb * (curl(dof(a)) - b)) * curl(tf(a)));
    
    // This leads to a classical Newton iteration where the unknown is the da increment:
    magnetostatics += integral(wholedomain, nu * curl(dof(a)) * curl(tf(a)));
    // The contribution below has tag number 1 so it can be assembled individually.
    // Dropping this line leads to a fixed point iteration.
    magnetostatics += integral(steel, dhdb * curl(dof(a)) * curl(tf(a)), 0, 1);

    // Current density source. Use ports to set a total current flow condition instead.
    double js = 6e4;
    magnetostatics += integral(conductor1, -array3x1(0,0,js) * tf(a));
    magnetostatics += integral(conductor2, -array3x1(0,0,-js) * tf(a));

    std::cout << "Total current in each conductor is " << js*expression(1).integrate(conductor1, 5) << " A" << std::endl;

    // Initial solution is a = 0 (thus b = 0):
    vec x(magnetostatics);
    
    // Newton iteration to solve the nonlinear problem A(x)*x = rhs(x).
    // At each iteration we solve J(x)*dx = rhs(x) - A(x)*x where J is the Jacobian matrix.
    double normdx = 1, maxb; int iter = 0;
    while (normdx > 1e-10)
    {
        // Generate the A and rhs terms (they all have tag 0 by default):
        magnetostatics.generate(0);

        // Get A and leave the generated terms in the formulation for reuse when getting J:
        mat A = magnetostatics.A(true);
        vec rhs = magnetostatics.b();
        
        // Generate the extra contribution of tag 1 to get J (J = A + sensitivity terms):
        magnetostatics.generate(1);
        // Get J and clear the generated terms:
        mat J = magnetostatics.A();

        vec dx = solve(J, rhs - A*x);
       
        x = x + dx;
        setdata(x);

        // The stopping criterion is a small enough increment:
        normdx = dx.norm();
        std::cout << "Increment norm @" << iter << " is " << normdx;
        
        // Do not allow b to go out of the provided [0, 2.4] T data range for the h(b) curve:
        maxb = norm(b).max(wholedomain, 5)[0];
        std::cout << " (max b is " << maxb << " T)" << std::endl;
        if (maxb > 2.4)
            x = 2.4/maxb * x;
        
        setdata(x);
        
        b.write(wholedomain, "b.vtu", 2);
        
        iter++;
    }
    
    // Code validation line. Can be removed.
    std::cout << (std::abs(maxb - 1.8239373092)/1.8239373092 < 1e-10);
}

