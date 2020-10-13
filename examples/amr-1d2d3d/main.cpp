// The purpose of this example is to check that AMR works as expected when the software is modified.


#include "sparselizardbase.h"


using namespace mathop;

double hadapt1d(void)
{
    int line = 1, left = 2, right = 3;

    mesh mymesh("1d.msh");

    field v("h1"), x("x");

    mymesh.setadaptivity(x, 0, 2);

    v.setorder(line, 1);

    v.setconstraint(left, 10);
    v.setconstraint(right, 2);

    formulation poisson;

    poisson += integral(line, grad(dof(v))*grad(tf(v)));

    solve(poisson);

    adapt(1);
    adapt(1);

    double intgr = norm(grad(v)).integrate(line, 5);
    double exact = 8.0;
    double relerror = std::abs(intgr-exact)/exact;

    return relerror;
}

double hadapt2d(void)
{	
    int s1 = 1, s2 = 2, s3 = 3, s4 = 4, l1 = 5, l2 = 6, l3 = 7, l4 = 8;

    mesh mymesh("2d.msh");

    int sur = regionunion({s1,s2,s3,s4});
    int lin = regionunion({l1,l2,l3,l4});

    field v("h1"), x("x"), y("y");

    expression criterion = 1+ifpositive(sin(3*x*y), sin(50*x)*sin(57*y)*x*y, 0);

    mymesh.setadaptivity(criterion, 0, 5);

    v.setorder(s1, 5);
    v.setorder(s2, 4);
    v.setorder(s3, 3);
    v.setorder(s4, 3);

    v.setconstraint(lin, x);

    formulation poisson;

    poisson += integral(sur, grad(dof(v))*grad(tf(v)));

    for (int i = 0; i < 5; i++)
        adapt(1);
    
    solve(poisson);

    double intgr = norm(grad(v)).integrate(sur, 5);
    double exact = 2.0*getpi();
    double relerror = std::abs(intgr-exact)/exact;

    return relerror;
}

double hadapt3d(void)
{	
    int vol = 1, left = 2, right = 3;

    mesh mymesh("3d.msh");

    field v("h1"), x("x"), y("y"), z("z");
    v.setorder(vol, 1);

    expression criterion = 1+ifpositive(sin(5*x)*sin(4*y)*sin(6*z), sin(50*x)*sin(57*y)*y*sin(53*z), 0);

    mymesh.setadaptivity(criterion, 0, 5);

    v.setconstraint(left, 10);
    v.setconstraint(right, 0);

    formulation poisson;

    poisson += integral(vol, grad(dof(v))*grad(tf(v)));

    for (int i = 0; i < 5; i++)
        adapt(1);
    
    solve(poisson);

    double intgr = norm(grad(v)).integrate(vol, 5);
    double exact = 10.0;
    double relerror = std::abs(intgr-exact)/exact;

    return relerror;
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    double relerror1d = hadapt1d();
    double relerror2d = hadapt2d();
    double relerror3d = hadapt3d();

    std::cout << relerror1d << " " << relerror2d << " " << relerror3d << std::endl;

    // Code validation line. Can be removed.
    std::cout << (relerror1d < 2e-15 && relerror2d < 4e-11 && relerror3d < 2e-12);

    SlepcFinalize();

    return 0;
}

