#include "opfieldorder.h"


opfieldorder::opfieldorder(std::vector<std::shared_ptr<rawfield>> fieldsin, double alpha, double absthres)
{
    myfields = fieldsin;
    myalpha = alpha;
    mythreshold = absthres;
}

std::vector<std::vector<densemat>> opfieldorder::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    // Make sure all fields have the same interpolation order:
    for (int i = 1; i < myfields.size(); i++)
    {
        if (myfields[i]->getinterpolationorders() != myfields[0]->getinterpolationorders())
        {
            logs log;
            log.msg() << "Error in 'opfieldorder' object: this object expects the same interpolation order for all subfields" << std::endl;
            log.error();
        }
    }
    
    int numelems = elemselect.countinselection();
    int elementtypenumber = elemselect.getelementtypenumber();
    std::vector<int> elems = elemselect.getelementnumbers();
    std::vector<int> fieldorders;
    int maxorder = myfields[0]->getinterpolationorders(elementtypenumber, elems, fieldorders);
    
    if (myalpha != -1.0)
    {
        indexmat fo(numelems, 1, fieldorders);
        std::vector<std::vector<int>> splitorders = fo.findalloccurences(maxorder);
        
        for (int o = 0; o <= maxorder; o++)
        {
            int numinorder = splitorders[o].size();
            if (numinorder == 0)
                continue;
                
            std::vector<int> elemsinorder(numinorder);
            for (int i = 0; i < numinorder; i++)
                elemsinorder[i] = elems[splitorders[o][i]];
        
            // Sum the weights of each field:
            std::vector<double> weightsforeachorder;
            myfields[0]->getweightsforeachorder(elementtypenumber, o, elemsinorder, weightsforeachorder);
            for (int i = 1; i < myfields.size(); i++)
            {
                std::vector<double> curweightsforeachorder;
                myfields[i]->getweightsforeachorder(elementtypenumber, o, elemsinorder, curweightsforeachorder);
                for (int j = 0; j < weightsforeachorder.size(); j++)
                    weightsforeachorder[j] += curweightsforeachorder[j];
            }
            
            std::vector<int> lowestorders;
            myfields[0]->getinterpolationorders(o, myalpha, mythreshold, weightsforeachorder, lowestorders);

            for (int i = 0; i < numinorder; i++)
                fieldorders[splitorders[o][i]] = lowestorders[i];
        }
    }
    
    densemat output(numelems, 1);
    double* outputvals = output.getvalues();
    for (int i = 0; i < numelems; i++)
        outputvals[i] = fieldorders[i];
    
    output = output.duplicatehorizontally(evaluationcoordinates.size()/3);
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), {{},{output}});

    // The field order is on the cos0 harmonic:
    return {{},{output}};
}

densemat opfieldorder::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }

    densemat output = interpolate(elemselect, evaluationcoordinates, meshdeform)[1][0];

    output = output.getflattened();
    output = output.duplicatevertically(numtimeevals);

    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
        
    return output;
}

std::shared_ptr<operation> opfieldorder::copy(void)
{
    std::shared_ptr<opfieldorder> op(new opfieldorder(myfields, myalpha, mythreshold));
    *op = *this;
    op->reuse = false;
    return op;
}

void opfieldorder::print(void)
{
    std::cout << "fieldorder";
}

