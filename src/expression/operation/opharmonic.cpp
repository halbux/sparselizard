#include "opharmonic.h"


std::vector<std::vector<densematrix>> opharmonic::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<std::vector<densematrix>> argmat;
    if (mynumfftharms < 0)
        argmat = myarg->interpolate(elemselect, evaluationcoordinates, meshdeform);
    else
    {
        densematrix timevals = myarg->multiharmonicinterpolate(mynumfftharms, elemselect, evaluationcoordinates, meshdeform);
        argmat = myfft::fft(timevals, elemselect.countinselection(), evaluationcoordinates.size()/3);
    }

    densematrix output;
    if (argmat.size() > myharmnum && argmat[myharmnum].size() > 0)
        output = argmat[myharmnum][0];
    else
        output = densematrix(elemselect.countinselection(), evaluationcoordinates.size()/3, 0);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), {{},{output}});
    
    return {{},{output}};
}

densematrix opharmonic::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densematrix output = interpolate(elemselect, evaluationcoordinates, meshdeform)[1][0];
    output = output.flatten();
    output = output.duplicatevertically(numtimeevals);
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opharmonic::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    return shared_from_this();
}

std::shared_ptr<operation> opharmonic::copy(void)
{
    std::shared_ptr<opharmonic> op(new opharmonic(myharmnum, myarg, mynumfftharms));
    *op = *this;
    op->reuse = false;
    return op;
}

void opharmonic::print(void)
{
    std::cout << "harmonic(" << myharmnum << ", ";
    myarg->print();
    std::cout << ")";
}
