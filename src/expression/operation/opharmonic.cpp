#include "opharmonic.h"


opharmonic::opharmonic(std::vector<int> origharms, std::vector<int> destharms, std::shared_ptr<operation> arg, int numfftharms)
{
    myorigharms = origharms;
    mydestharms = destharms;
    myarg = arg;
    mynumfftharms = numfftharms;
}
    
std::vector<std::vector<densematrix>> opharmonic::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    int numelems = elemselect.countinselection();
    
    std::vector<std::vector<densematrix>> argmat;
    if (mynumfftharms < 0)
        argmat = myarg->interpolate(elemselect, evaluationcoordinates, meshdeform);
    else
    {
        densematrix timevals = myarg->multiharmonicinterpolate(mynumfftharms, elemselect, evaluationcoordinates, meshdeform);
        argmat = myfft::fft(timevals, numelems, evaluationcoordinates.size()/3);
    }
    
    int maxdestharm = *std::max_element(mydestharms.begin(), mydestharms.end());
    std::vector<std::vector<densematrix>> output(maxdestharm+1, std::vector<densematrix>(0));

    for (int i = 0; i < myorigharms.size(); i++)
    {
        int horig = myorigharms[i];
        int hdest = mydestharms[i];
  
        if (horig < argmat.size() && argmat[horig].size() > 0)
        {
            if (output[hdest].size() == 0)
                output[hdest] = {argmat[horig][0]};
            else
                output[hdest][0].add(argmat[horig][0]);
        }
        else
        {
            if (output[hdest].size() == 0)
                output[hdest] = {densematrix(numelems, evaluationcoordinates.size()/3, 0)};
        }
    }

    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), output);
    
    return output;
}

densematrix opharmonic::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    std::vector<std::vector<densematrix>> interpolated = interpolate(elemselect, evaluationcoordinates, meshdeform);
    // Compute at 'numtimevals' instants in time the multiharmonic data:
    densematrix output = myfft::inversefft(interpolated, numtimeevals, elemselect.countinselection(), evaluationcoordinates.size()/3);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

bool opharmonic::isharmonicone(std::vector<int> disjregs)
{
    for (int i = 0; i < mydestharms.size(); i++)
    {
        if (mydestharms[i] != 1)
            return false;
    }
    return true;
}

std::shared_ptr<operation> opharmonic::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    return shared_from_this();
}

std::shared_ptr<operation> opharmonic::copy(void)
{
    std::shared_ptr<opharmonic> op(new opharmonic(myorigharms, mydestharms, myarg, mynumfftharms));
    *op = *this;
    op->reuse = false;
    return op;
}

void opharmonic::print(void)
{
    std::cout << "moveharmonic(";
    myarg->print();
    std::cout << ")";
}
