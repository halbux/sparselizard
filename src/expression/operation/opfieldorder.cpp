#include "opfieldorder.h"


std::vector<std::vector<densematrix>> opfieldorder::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<int> drsinelsel = elemselect.getdisjointregions();
    
    densematrix output(elemselect.countinselection(), 1);
    
    elements* myelems = universe::mymesh->getelements();
    int numdisjregs = universe::mymesh->getdisjointregions()->count();
    std::vector<int> orderofdisjregs(numdisjregs, -1);
    for (int i = 0; i < drsinelsel.size(); i++)
        orderofdisjregs[drsinelsel[i]] = myfield->getinterpolationorder(drsinelsel[i]);
    
    int elementtypenumber = elemselect.getelementtypenumber();
    
    std::vector<int> elems = elemselect.getelementnumbers();
    
    double* vals = output.getvalues();
    for (int i = 0; i < output.countrows(); i++)
        vals[i] = orderofdisjregs[ myelems->getdisjointregion(elementtypenumber, elems[i]) ];
    
    output = output.duplicatehorizontally(evaluationcoordinates.size()/3);
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), {{},{output}});

    // The field order is on the cos0 harmonic:
    return {{},{output}};
}

densematrix opfieldorder::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    std::vector<int> drsinelsel = elemselect.getdisjointregions();
    
    densematrix output(elemselect.countinselection(), 1);
    
    elements* myelems = universe::mymesh->getelements();
    int numdisjregs = universe::mymesh->getdisjointregions()->count();
    std::vector<int> orderofdisjregs(numdisjregs, -1);
    for (int i = 0; i < drsinelsel.size(); i++)
        orderofdisjregs[drsinelsel[i]] = myfield->getinterpolationorder(drsinelsel[i]);
    
    int elementtypenumber = elemselect.getelementtypenumber();
    
    std::vector<int> elems = elemselect.getelementnumbers();
    
    double* vals = output.getvalues();
    for (int i = 0; i < output.countrows(); i++)
        vals[i] = orderofdisjregs[ myelems->getdisjointregion(elementtypenumber, elems[i]) ];
    
    output = output.duplicatehorizontally(evaluationcoordinates.size()/3);
    output = output.flatten();
    output = output.duplicatevertically(numtimeevals);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
        
    return output;
}

std::shared_ptr<operation> opfieldorder::copy(void)
{
    std::shared_ptr<opfieldorder> op(new opfieldorder(myfield));
    *op = *this;
    op->reuse = false;
    return op;
}

void opfieldorder::print(void) { std::cout << "fieldorder"; }

