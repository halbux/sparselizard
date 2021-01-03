#include "opdetjac.h"


std::vector<std::vector<densematrix>> opdetjac::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Compute the Jacobian terms or reuse if available in the universe.
    std::shared_ptr<jacobian> myjac;
    if (universe::isreuseallowed && universe::computedjacobian != NULL)
        myjac = universe::computedjacobian;
    else
        myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));
    
    if (universe::isreuseallowed)
        universe::computedjacobian = myjac;
    
    // The detjac is on the cos0 harmonic:
    return {{},{myjac->getdetjac().copy()}};
}

densematrix opdetjac::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Compute the Jacobian terms or reuse if available in the universe.
    std::shared_ptr<jacobian> myjac;
    if (universe::isreuseallowed && universe::computedjacobian != NULL)
        myjac = universe::computedjacobian;
    else
        myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));
    
    if (universe::isreuseallowed)
        universe::computedjacobian = myjac;
    
    densematrix computeddetjac = (myjac->getdetjac().copy());
    
    computeddetjac = computeddetjac.getflattened();
    return computeddetjac.duplicatevertically(numtimeevals);
}

std::shared_ptr<operation> opdetjac::copy(void)
{
    std::shared_ptr<opdetjac> op(new opdetjac);
    *op = *this;
    return op;
}

void opdetjac::print(void) { std::cout << "detjac"; }



