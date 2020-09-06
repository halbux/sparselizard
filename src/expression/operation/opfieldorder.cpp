#include "opfieldorder.h"


opfieldorder::opfieldorder(std::shared_ptr<rawfield> fieldin, double alpha, double absthres)
{
    myfield = fieldin;
    myalpha = alpha;
    mythreshold = absthres;
}

std::vector<std::vector<densematrix>> opfieldorder::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    int numelems = elemselect.countinselection();
    int elementtypenumber = elemselect.getelementtypenumber();
    std::vector<int> elems = elemselect.getelementnumbers();
    std::vector<int> fieldorders;
    int maxorder = myfield->getinterpolationorders(elementtypenumber, elems, fieldorders);
    
    if (myalpha != -1.0)
    {
        intdensematrix fo(numelems, 1, fieldorders);
        std::vector<std::vector<int>> splitorders = fo.findalloccurences(maxorder);
        
        for (int o = 0; o <= maxorder; o++)
        {
            int numinorder = splitorders[o].size();
            if (numinorder == 0)
                continue;
                
            std::vector<int> elemsinorder(numinorder);
            for (int i = 0; i < numinorder; i++)
                elemsinorder[i] = elems[splitorders[o][i]];
        
            std::vector<int> lowestorders;
            myfield->getinterpolationorders(elementtypenumber, o, elemsinorder, myalpha, mythreshold, lowestorders);

            for (int i = 0; i < numinorder; i++)
                fieldorders[splitorders[o][i]] = lowestorders[i];
        }
    }
    
    densematrix output(numelems, 1);
    double* outputvals = output.getvalues();
    for (int i = 0; i < numelems; i++)
        outputvals[i] = fieldorders[i];
    
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

    int numelems = elemselect.countinselection();
    int elementtypenumber = elemselect.getelementtypenumber();
    std::vector<int> elems = elemselect.getelementnumbers();
    std::vector<int> fieldorders;
    int maxorder = myfield->getinterpolationorders(elementtypenumber, elems, fieldorders);
    
    if (myalpha != -1.0)
    {
        intdensematrix fo(numelems, 1, fieldorders);
        std::vector<std::vector<int>> splitorders = fo.findalloccurences(maxorder);
        
        for (int o = 0; o <= maxorder; o++)
        {
            int numinorder = splitorders[o].size();
            if (numinorder == 0)
                continue;
                
            std::vector<int> elemsinorder(numinorder);
            for (int i = 0; i < numinorder; i++)
                elemsinorder[i] = elems[splitorders[o][i]];
        
            std::vector<int> lowestorders;
            myfield->getinterpolationorders(elementtypenumber, o, elemsinorder, myalpha, mythreshold, lowestorders);

            for (int i = 0; i < numinorder; i++)
                fieldorders[splitorders[o][i]] = lowestorders[i];
        }
    }
    
    densematrix output(numelems, 1);
    double* outputvals = output.getvalues();
    for (int i = 0; i < numelems; i++)
        outputvals[i] = fieldorders[i];
    
    output = output.duplicatehorizontally(evaluationcoordinates.size()/3);
    output = output.flatten();
    output = output.duplicatevertically(numtimeevals);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
        
    return output;
}

std::shared_ptr<operation> opfieldorder::copy(void)
{
    std::shared_ptr<opfieldorder> op(new opfieldorder(myfield, myalpha, mythreshold));
    *op = *this;
    op->reuse = false;
    return op;
}

void opfieldorder::print(void)
{
    std::cout << "fieldorder(";
    myfield->print();
    std::cout << ")";
}

