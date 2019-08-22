#include "opparameter.h"


std::vector<std::vector<densematrix>> opparameter::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<std::vector<densematrix>> output = myparameter->interpolate(myrow, mycolumn, elemselect, evaluationcoordinates, meshdeform);
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), output);
    return output;
}

densematrix opparameter::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densematrix output = myparameter->multiharmonicinterpolate(myrow, mycolumn, numtimeevals, elemselect, evaluationcoordinates, meshdeform);
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    return output;
}

bool opparameter::isharmonicone(std::vector<int> disjregs) 
{ 
    for (int i = 0; i < disjregs.size(); i++)
    {
        if ( (myparameter->get(disjregs[i], myrow, mycolumn))->isharmonicone({disjregs[i]}) == false )
            return false;
    }
    return true;        
}

std::shared_ptr<operation> opparameter::simplify(std::vector<int> disjregs)
{
    for (int i = 0; i < disjregs.size(); i++)
        myparameter->simplify(myrow, mycolumn, disjregs[i]);
    return shared_from_this();
}

bool opparameter::isvalueorientationdependent(std::vector<int> disjregs)
{
    for (int i = 0; i < disjregs.size(); i++)
    {
        if ( (myparameter->get(disjregs[i], myrow, mycolumn))->isvalueorientationdependent({disjregs[i]}) == true )
            return true;
    }
    return false;
}

std::shared_ptr<operation> opparameter::copy(void)
{
    std::shared_ptr<opparameter> op(new opparameter(myparameter, myrow, mycolumn));
    *op = *this;
    return op;
}

void opparameter::print(void) { std::cout << "param" << myparameter->countrows() << "x" << myparameter->countcolumns() << "(" << myrow << "," << mycolumn << ")"; }
