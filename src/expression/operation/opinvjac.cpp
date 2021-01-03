#include "opinvjac.h"


std::vector<std::vector<densematrix>> opinvjac::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Compute the Jacobian terms or reuse if available in the universe.
    std::shared_ptr<jacobian> myjac;
    if (universe::isreuseallowed && universe::computedjacobian != NULL)
        myjac = universe::computedjacobian;
    else
        myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));
    
    if (universe::isreuseallowed)
        universe::computedjacobian = myjac;
    
    // The invjac is on the cos0 harmonic:
    return {{},{(myjac->getinvjac(myrow,mycol)).copy()}};
}

densematrix opinvjac::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Compute the Jacobian terms or reuse if available in the universe.
    std::shared_ptr<jacobian> myjac;
    if (universe::isreuseallowed && universe::computedjacobian != NULL)
        myjac = universe::computedjacobian;
    else
        myjac = std::shared_ptr<jacobian>(new jacobian(elemselect, evaluationcoordinates, meshdeform));
    
    if (universe::isreuseallowed)
        universe::computedjacobian = myjac;
    
    densematrix computedinvjac = (myjac->getinvjac(myrow,mycol)).copy();
    
    computedinvjac = computedinvjac.getflattened();
    return computedinvjac.duplicatevertically(numtimeevals);
}

std::shared_ptr<operation> opinvjac::copy(void)
{
    std::shared_ptr<opinvjac> op(new opinvjac(myrow, mycol));
    *op = *this;
    return op;
}

void opinvjac::print(void) { std::cout << "invjac" << myrow << mycol; }



