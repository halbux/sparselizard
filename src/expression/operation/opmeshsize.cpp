#include "opmeshsize.h"


std::vector<std::vector<densematrix>> opmeshsize::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    int typenum = elemselect.getelementtypenumber();
    gausspoints mygp(typenum, myintegrationorder);

    int numgp = mygp.count();
    std::vector<double> evalcoords = mygp.getcoordinates();
    std::vector<double> weights = mygp.getweights();
    densematrix weightsmat(numgp,1, weights);

    std::shared_ptr<opdetjac> op(new opdetjac);

    bool wasreuseallowed = universe::isreuseallowed;
    universe::isreuseallowed = false;
    densematrix output = op->interpolate(elemselect, evalcoords, meshdeform)[1][0];
    universe::isreuseallowed = wasreuseallowed;
    
    output.abs();
    output = output.multiply(weightsmat);
    output = output.duplicatehorizontally(evaluationcoordinates.size()/3);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), {{},{output}});

    // The mesh size is on the cos0 harmonic:
    return {{},{output}};
}

densematrix opmeshsize::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    int typenum = elemselect.getelementtypenumber();
    gausspoints mygp(typenum, myintegrationorder);

    int numgp = mygp.count();
    std::vector<double> evalcoords = mygp.getcoordinates();
    std::vector<double> weights = mygp.getweights();
    densematrix weightsmat(numgp,1, weights);

    std::shared_ptr<opdetjac> op(new opdetjac);

    bool wasreuseallowed = universe::isreuseallowed;
    universe::isreuseallowed = false;
    densematrix output = op->interpolate(elemselect, evalcoords, meshdeform)[1][0];
    universe::isreuseallowed = wasreuseallowed;
    
    output.abs();
    output = output.multiply(weightsmat);
    output = output.duplicatehorizontally(evaluationcoordinates.size()/3);
    output = output.getflattened();
    output = output.duplicatevertically(numtimeevals);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
        
    return output;
}

std::shared_ptr<operation> opmeshsize::copy(void)
{
    std::shared_ptr<opmeshsize> op(new opmeshsize(myintegrationorder));
    *op = *this;
    op->reuse = false;
    return op;
}

void opmeshsize::print(void) { std::cout << "meshsize"; }

