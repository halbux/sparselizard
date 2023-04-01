#include "opabs.h"


std::vector<std::vector<densemat>> opabs::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<std::vector<densemat>> argmat = myarg->interpolate(elemselect, evaluationcoordinates, meshdeform);
    
    if (argmat.size() == 2 && argmat[1].size() == 1)
    {
        argmat[1][0].abs();
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), argmat);
        
        return argmat;
    }

    logs log;
    log.msg() << "Error in 'opabs' object: without FFT abs() can only be computed for constant (harmonic 1) operations" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat opabs::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat output = myarg->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    output.abs();
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opabs::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    if (myarg->isconstant())
        return std::shared_ptr<operation>(new opconstant(std::abs(myarg->getvalue())));
    else
        return shared_from_this();
}

std::shared_ptr<operation> opabs::copy(void)
{
    std::shared_ptr<opabs> op(new opabs(myarg));
    *op = *this;
    op->reuse = false;
    return op;
}

double opabs::evaluate(void)
{
    return std::abs(myarg->evaluate());
}

std::vector<double> opabs::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluated = myarg->evaluate(xcoords, ycoords, zcoords);
    for (int i = 0; i < evaluated.size(); i++)
        evaluated[i] = std::abs(evaluated[i]);
    return evaluated;
}

void opabs::print(void)
{
    std::cout << "abs(";
    myarg->print();
    std::cout << ")";
}
