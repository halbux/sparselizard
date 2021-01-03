#include "opjac.h"


std::vector<std::vector<densematrix>> opjac::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Compute the Jacobian terms or reuse if available in the universe.
    std::shared_ptr<jacobian> myjac;
    if (universe::isreuseallowed && universe::computedjacobian != NULL)
        myjac = universe::computedjacobian;
    else
        myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));
    
    if (universe::isreuseallowed)
        universe::computedjacobian = myjac;
    
    // The jac is on the cos0 harmonic:
    return {{},{(myjac->getjac(myrow,mycol)).copy()}};
}

densematrix opjac::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Compute the Jacobian terms or reuse if available in the universe.
    std::shared_ptr<jacobian> myjac;
    if (universe::isreuseallowed && universe::computedjacobian != NULL)
        myjac = universe::computedjacobian;
    else
        myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));
    
    if (universe::isreuseallowed)
        universe::computedjacobian = myjac;
    
    densematrix computedjac = (myjac->getjac(myrow,mycol)).copy();
    
    computedjac = computedjac.getflattened();
    return computedjac.duplicatevertically(numtimeevals);
}

std::shared_ptr<operation> opjac::copy(void)
{
    std::shared_ptr<opjac> op(new opjac(myrow, mycol));
    *op = *this;
    return op;
}

void opjac::print(void) { std::cout << "jac" << myrow << mycol; }



