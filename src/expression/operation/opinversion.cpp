#include "opinversion.h"


std::vector<std::vector<densematrix>> opinversion::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
	std::vector<std::vector<densematrix>> argmat = myarg->interpolate(elemselect, evaluationcoordinates, meshdeform);
	
    if (argmat.size() == 2 && argmat[1].size() == 1)
    {
        argmat[1][0].invert();
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), argmat);
        
        return argmat;
    }

    std::cout << "Error in 'opinversion' object: without FFT divisions can only be computed for constant (harmonic 1) operations" << std::endl;
    abort();
}

densematrix opinversion::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densematrix output = myarg->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    output.invert();
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opinversion::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    if (myarg->isconstant())
        return std::shared_ptr<operation>(new opconstant(1.0/myarg->getvalue()));
    else
        return shared_from_this();
}

std::shared_ptr<operation> opinversion::copy(void)
{
    std::shared_ptr<opinversion> op(new opinversion(myarg));
    *op = *this;
    op->reuse = false;
    return op;
}

void opinversion::print(void)
{
    std::cout << "1/(";
    myarg->print();
    std::cout << ")";
}
