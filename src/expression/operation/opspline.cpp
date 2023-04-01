#include "opspline.h"


opspline::opspline(spline spl, std::shared_ptr<operation> arg)
{
    myarg = arg; myspline = spl; 
}
        
std::vector<std::vector<densemat>> opspline::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
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
        argmat[1][0] = myspline.evalat(argmat[1][0]);
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), argmat);
        
        return argmat;
    }

    logs log;
    log.msg() << "Error in 'opspline' object: without FFT a spline can only be interpolated for constant (harmonic 1) operations" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat opspline::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat output = myarg->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    output = myspline.evalat(output);
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opspline::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    return shared_from_this();
}

std::shared_ptr<operation> opspline::copy(void)
{
    std::shared_ptr<opspline> op(new opspline(myspline, myarg));
    *op = *this;
    op->reuse = false;
    return op;
}

void opspline::print(void)
{
    std::cout << "spline(";
    myarg->print();
    std::cout << ")";
}
