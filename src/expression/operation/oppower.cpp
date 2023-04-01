#include "oppower.h"


std::vector<std::vector<densemat>> oppower::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<std::vector<densemat>> computedbase = mybase->interpolate(elemselect, evaluationcoordinates, meshdeform);
    std::vector<std::vector<densemat>> computedexponent = myexponent->interpolate(elemselect, evaluationcoordinates, meshdeform);

    if (computedbase.size() == 2 && computedbase[1].size() == 1 && computedexponent.size() == 2 && computedexponent[1].size() == 1)
    {
        computedbase[1][0].power(computedexponent[1][0]);
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), computedbase);
        
        return computedbase;
    }

    logs log;
    log.msg() << "Error in 'oppower' object: without FFT a power can only be computed for constant (harmonic 1) operations. Use an FFT call or rewrite as a product if possible." << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat oppower::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat computedbase = mybase->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    densemat computedexponent = myexponent->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);

    computedbase.power(computedexponent);
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), computedbase);
    
    return computedbase;
}

std::shared_ptr<operation> oppower::simplify(std::vector<int> disjregs)
{
    mybase = mybase->simplify(disjregs);
    myexponent = myexponent->simplify(disjregs);

    if (mybase->isconstant() && myexponent->isconstant())
        return std::shared_ptr<operation>(new opconstant(std::pow(mybase->getvalue(),myexponent->getvalue())));
    else
        return shared_from_this();
}

std::shared_ptr<operation> oppower::copy(void)
{
    std::shared_ptr<oppower> op(new oppower(mybase, myexponent));
    *op = *this;
    op->reuse = false;
    return op;
}

double oppower::evaluate(void)
{
    return std::pow(mybase->evaluate(), myexponent->evaluate());
}

std::vector<double> oppower::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluatedbase = mybase->evaluate(xcoords, ycoords, zcoords);
    std::vector<double> evaluatedexponent = myexponent->evaluate(xcoords, ycoords, zcoords);

    for (int i = 0; i < xcoords.size(); i++)
        evaluatedbase[i] = std::pow(evaluatedbase[i], evaluatedexponent[i]);
    return evaluatedbase;
}

void oppower::print(void)
{
    std::cout << "(";
    mybase->print();
    std::cout << ")^(";
    myexponent->print();
    std::cout << ")";
}
