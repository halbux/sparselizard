#include "opsum.h"


void opsum::subtractterm(std::shared_ptr<operation> term)
{
    std::shared_ptr<opproduct> opprod(new opproduct);
    std::shared_ptr<opconstant> opconst(new opconstant(-1));
    
    opprod->multiplybyterm(term);
    opprod->multiplybyterm(opconst);
    
    sumterms.push_back(opprod);
}

std::vector<std::vector<densemat>> opsum::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    // Compute first all the sum terms and get the max number of harmonics.
    int maxnumberofharmonics = 0;
    std::vector< std::vector<std::vector<densemat>> > computedterms(sumterms.size());
    for (int i = 0; i < sumterms.size(); i++)
    {
        computedterms[i] = sumterms[i]->interpolate(elemselect, evaluationcoordinates, meshdeform);
        if (computedterms[i].size() > maxnumberofharmonics)
            maxnumberofharmonics = computedterms[i].size();
    }
    
    // Initialise the output vector holding the sum:
    std::vector<std::vector<densemat>> output(maxnumberofharmonics, std::vector<densemat> {});
    
    // Sum all harmonics together (ignore harmonic sin0):
    for (int harm = 1; harm < maxnumberofharmonics; harm++)
    {
        for (int i = 0; i < sumterms.size(); i++)
        {
            // If the sum term exists for the current harmonic:
            if (harm < computedterms[i].size() && computedterms[i][harm].size() == 1)
            {
                // If the term does not exist in the output:
                if (output[harm].size() == 0)
                    output[harm] = {computedterms[i][harm][0]};
                else
                    output[harm][0].add(computedterms[i][harm][0]);
            }
        }
    }
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), output);
    
    return output;
}

densemat opsum::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat output = sumterms[0]->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    
    for (int i = 1; i < sumterms.size(); i++)
        output.add(sumterms[i]->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform));
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opsum::expand(void)
{
    for (int i = 0; i < sumterms.size(); i++)
    {
        // We only want to expand operations that include a dof(), tf() or port.
        if (sumterms[i]->isdofincluded() || sumterms[i]->istfincluded() || sumterms[i]->isportincluded())
            sumterms[i] = sumterms[i]->expand();
    }
    // Regroup all sum terms in this one:
    group();
    
    return shared_from_this();
}

void opsum::group(void)
{
    std::vector<std::shared_ptr<operation>> groupedsumterms = {};
    
    for (int i = 0; i < sumterms.size(); i++)
    {
        sumterms[i]->group();

        // We want to keep the operations to reuse untouched!
        if (sumterms[i]->issum() && sumterms[i]->isreused() == false)
        {
            for (int j = 0; j < sumterms[i]->count(); j++)
                groupedsumterms.push_back(sumterms[i]->getargument(j));
        }
        else
            groupedsumterms.push_back(sumterms[i]);
    }
    sumterms = groupedsumterms;
}

std::shared_ptr<operation> opsum::simplify(std::vector<int> disjregs)
{
    for (int i = 0; i < count(); i++)
        sumterms[i] = sumterms[i]->simplify(disjregs);
    
    group();

    double summed = 0;
    int numsumterms = count();
    for (int i = numsumterms-1; i >= 0; i--)
    {
        if (sumterms[i]->isconstant())
        {
            summed += sumterms[i]->getvalue();
            removeterm(i);
        }
    }
    // If there were not only constants in the sum:
    if (count() > 0)
    {
        if (summed != 0)
            addterm(std::shared_ptr<operation>(new opconstant(summed)));
    }
    else
    {
        // We do not want to invalidate this operation and thus add the constant term:
        sumterms = {std::shared_ptr<operation>(new opconstant(summed))};
        return sumterms[0];
    }
    
    return shared_from_this();
}

std::shared_ptr<operation> opsum::copy(void)
{
    std::shared_ptr<opsum> op(new opsum);
    *op = *this;
    op->reuse = false;
    return op;
}

std::vector<double> opsum::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluated(xcoords.size(), 0);
    for (int i = 0; i < sumterms.size(); i++)
    {
        std::vector<double> current = sumterms[i]->evaluate(xcoords, ycoords, zcoords);
        for (int j = 0; j < xcoords.size(); j++)
            evaluated[j] += current[j];
    }
    return evaluated;
}

void opsum::print(void)
{
    std::cout << "(";
    for (int i = 0; i < sumterms.size(); i++)
    {
        if (i > 0)
            std::cout << " + ";
        sumterms[i]->print();
    }
    std::cout << ")";
}
