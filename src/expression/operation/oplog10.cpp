#include "oplog10.h"


std::vector<std::vector<densematrix>> oplog10::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
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
        argmat[1][0].log10();
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), argmat);
        
        return argmat;
    }

    std::cout << "Error in 'oplog10' object: without FFT log10() can only be computed for constant (harmonic 1) operations" << std::endl;
    abort();
}

densematrix oplog10::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densematrix output = myarg->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    output.log10();

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> oplog10::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    if (myarg->isconstant())
        return std::shared_ptr<operation>(new opconstant(std::log10(myarg->getvalue())));
    else
        return shared_from_this();
}

std::shared_ptr<operation> oplog10::copy(void)
{
    std::shared_ptr<oplog10> op(new oplog10(myarg));
    *op = *this;
    op->reuse = false;
    return op;
}

double oplog10::evaluate(void)
{
    return std::log10(myarg->evaluate());
}

std::vector<double> oplog10::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluated = myarg->evaluate(xcoords, ycoords, zcoords);
    for (int i = 0; i < evaluated.size(); i++)
        evaluated[i] = std::log10(evaluated[i]);
    return evaluated;
}

void oplog10::print(void)
{
    std::cout << "log10(";
    myarg->print();
    std::cout << ")";
}
