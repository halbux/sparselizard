#include "opcos.h"


std::vector<std::vector<densemat>> opcos::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
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
        argmat[1][0].cos();
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), argmat);
        
        return argmat;
    }

    logs log;
    log.msg() << "Error in 'opcos' object: without FFT cos() can only be computed for constant (harmonic 1) operations" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat opcos::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat output = myarg->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    output.cos();
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opcos::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    if (myarg->isconstant())
        return std::shared_ptr<operation>(new opconstant(std::cos(myarg->getvalue())));
    else
        return shared_from_this();
}

std::shared_ptr<operation> opcos::copy(void)
{
    std::shared_ptr<opcos> op(new opcos(myarg));
    *op = *this;
    op->reuse = false;
    return op;
}

double opcos::evaluate(void)
{
    return std::cos(myarg->evaluate());
}

std::vector<double> opcos::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluated = myarg->evaluate(xcoords, ycoords, zcoords);
    for (int i = 0; i < evaluated.size(); i++)
        evaluated[i] = std::cos(evaluated[i]);
    return evaluated;
}

void opcos::print(void)
{
    std::cout << "cos(";
    myarg->print();
    std::cout << ")";
}
