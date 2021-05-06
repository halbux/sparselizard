// This is a benchmark example to validate 2D axisymmetry with 'hcurl' fields for electromagnetic waves.

#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    // The domain regions:
    int sur = 1, right = 2, boundary = 3;

    setaxisymmetry();

    int n = 50;

    shape q("quadrangle", sur, {0.005,0,0, 0.055,0,0, 0.055,0.04,0, 0.005,0.04,0}, {n,n,n,n});
    shape lb = q.getsons()[0];
    lb.setphysicalregion(boundary);
    shape lr = q.getsons()[1];
    lr.setphysicalregion(right);
    shape lt = q.getsons()[2];
    lt.setphysicalregion(boundary);
    shape ll = q.getsons()[3];
    ll.setphysicalregion(boundary);

    mesh mymesh({q,lb,lr,lt,ll});

    // Edge shape functions 'hcurl' for the electric field E.
    // Fields x and y are the x and y coordinate fields.
    field E("hcurl"), x("x"), y("y");

    // Use interpolation order 2 on the whole domain:
    E.setorder(sur, 2);

    double freq = 10.0e9, c = 299792458, k = 2.0*getpi()*freq/c;

    E.setconstraint(right, array3x1(0,y,0));
    // Perfect conductor boundary conditions:
    E.setconstraint(boundary);

    formulation maxwell;

    // This is the weak formulation for electromagnetic waves:
    maxwell += integral(sur, -curl(dof(E))*curl(tf(E)) + k*k*dof(E)*tf(E));

    maxwell.solve();  

    E.write(sur, "E.pos", 2);
    compx(E).write(sur, "Ex.pos", 2);
    compy(E).write(sur, "Ey.pos", 2);

    double maxE = norm(E).max(sur, 5)[0];

    std::cout << "Max electric field norm is " << maxE << " V/m" << std::endl;

    // Code validation line. Can be removed.
    std::cout << (maxE < 0.138395 && maxE > 0.138393);
}

