// The purpose of this example is to check that hp-adaptivity works as expected when the software is modified.


#include "sparselizardbase.h"


using namespace mathop;

double adapth1(void)
{	
    int sur = 1;
    
    int n = 5;
    shape q("quadrangle", sur, {0,0,0, 1,0,0, 1,1,0, 0,1,0}, {n,n,n,n});
    
    mesh mymesh({q}, 0);
    
    field u("h1"), v("h1"), w("h1"), x("x"), y("y");
    
    u.setorder(sur, 1);
    v.setorder(sur, 3);
    w.setorder(sur, 2);
    
    w.setvalue(sur, x*x);
    
    // P-adaptivity of field u and v:
    u.setorder(x, 2, 2);
    v.setorder(sin(5*x)*sin(3*y)+1e-3*u, 3, 5);
    
    adapt(0);
    adapt(0);
    
    u.setvalue(sur, x*x+y*y);
    
    // Increase the integration order during the hp-adaptivity field value projection:
    v.setupdateaccuracy(2);
    
    // H-adaptivity:
    mymesh.setadaptivity(sin(7*x)*y, 0, 3);
    
    // This can be exactly represented with an order 3 interpolation on triangles/quadrangles:
    expression toproject = x*x*x+y*y*y+0.25*y*x*x;
    
    formulation projection;

    projection += integral(sur, dof(v)*tf(v) - toproject*tf(v), 1);
    
    projection.generate();
    
    vec sol = solve(projection.A(), projection.b());
    
    // Adapt:
    while (adapt(0)) {};
    
    v.setdata(sur, sol);
    
    // Relative error on 'u' adaptation:
    double erroru = (abs(u).integrate(sur, 5)-2.0/3.0)/(2.0/3.0);
    // Relative error on 'sol' adaptation:
    double errorsol = abs(toproject-v).integrate(sur, 5)/abs(toproject).integrate(sur, 5);
    // Relative error on 'w' adaptation:
    double errorw = (abs(w).integrate(sur, 5)-1.0/3.0)/(1.0/3.0);
    
    erroru = std::abs(erroru); errorsol = std::abs(errorsol); errorw = std::abs(errorw);
    
    std::cout << "Relative error for 'h1': " << erroru << " " << errorsol << " " << errorw << std::endl;
    
    return erroru+errorsol+errorw;
}

double adapthcurl(void)
{	
    int vol = 1;
    
    mesh mymesh("box.msh", 0);
    
    field v("hcurl"), x("x"), y("y"), z("z");
    
    v.setorder(vol, 3);
    
    // This can be exactly represented with an order 3 interpolation on tetrahedra:
    expression toproject = array3x1(x*x*x, y*y*y, z*z*z);
    
    v.setvalue(vol, toproject, 1);
    
    // Increase the integration order during the h-adaptivity field value projection:
    v.setupdateaccuracy(1);
    
    // H-adaptivity:
    mymesh.setadaptivity(x*y, 0, 1);
    
    // Adapt:
    adapt(0);
    
    double errorv = norm(v-toproject).integrate(vol, 5)/norm(toproject).integrate(vol, 5);
    
    std::cout << "Relative error for 'hcurl': " << errorv << std::endl;
    
    return std::abs(errorv);
}

int main(void)
{	
    SlepcInitialize(0,{},0,0);

    double errorh1 = adapth1();
    double errorhcurl = adapthcurl();
    
    std::cout << (errorh1 < 3e-13 && errorhcurl < 2e-11);

    SlepcFinalize();

    return 0;
}

