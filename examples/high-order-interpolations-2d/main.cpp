// This code projects an expression on a finite element basis, i.e. it computes 
// the best approximation of the expression for the requested interpolation order.
// The relative L2 norm of the error is output so that one can investigate 
// the impact of the interpolation order and the mesh size on the error.
//
// The theory of finite elements states that for a fine enough mesh 
// (i.e. a large enough 'n' in the quad.geo file) the error should decrease 
// at a rate n^(p+1) where p is the interpolation order. 
// All standard shape functions in sparselizard follow this asymptotic decrease rate.


#include "sparselizard.h"


using namespace sl;

int main(void)
{	
    
    // CHANGE THE INTERPOLATION ORDER HERE TO SEE THE IMPACT ON THE ERROR
    // Recommended order range is 1 to 6 here. Note that when 
    // you reach machine precision the error stops decreasing.
    int interpolationorder = 5;
	

    // The domain region as defined in 'quad.geo':
    int face = 1;
    
    mesh mymesh("quad.msh");
    
    // Nodal shape functions 'h1' for field v.
    // Fields x and y are the x and y coordinate fields:
    field v("h1"), x("x"), y("y");

    // Use the requested interpolation order on the whole domain:
    v.setorder(face, interpolationorder);
    

    // A finite element approximation of this expression is computed:
    expression toproject = sin(10*x)*sin(13*y);
    
    // This formulation computes integral( v*v' ) = integral( toproject*v' ) : 
    formulation projection;
    
    // A last int argument is added. The 20 means the term below will be
    // assembled using an integration order +20 compared to the default
    // order (which is order of the dof + order of the tf + 2). 
    projection += integral(face, dof(v)*tf(v) - toproject*tf(v), +20 );
    
    projection.generate();

    vec solv = solve(projection.A(), projection.b());

    // Transfer the data from the solution vector to the v field:
    v.setdata(face, solv);
    
    // Compute the relative L2 norm of the error, i.e. the square root of 
    // the integral of the error squared divided by the expression squared.
    // The integration is exact for up to order 20 polynomials.
    double relativel2normoferror = sqrt( pow( v - toproject , 2).integrate(face, 20) / pow( toproject , 2).integrate(face, 20) );
    
    std::cout << std::endl << "Relative L2 norm of the error is " << relativel2normoferror << std::endl << std::endl;
    
    // Write with an order 4 interpolation:
    v.write(face, "v.pos", 4);
    
    // Code validation line. Can be removed.
    std::cout << (relativel2normoferror < 1e-5);
}

