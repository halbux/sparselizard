#include "opcondition.h"

opcondition::opcondition(std::shared_ptr<operation> condarg, std::shared_ptr<operation> truearg, std::shared_ptr<operation> falsearg)
{
    mycond = condarg;
    mytrue = truearg;
    myfalse = falsearg;
}


std::vector<std::vector<densemat>> opcondition::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<std::vector<densemat>> condargmat = mycond->interpolate(elemselect, evaluationcoordinates, meshdeform);
    std::vector<std::vector<densemat>> trueargmat = mytrue->interpolate(elemselect, evaluationcoordinates, meshdeform);
    std::vector<std::vector<densemat>> falseargmat = myfalse->interpolate(elemselect, evaluationcoordinates, meshdeform);

    if (condargmat.size() == 2 && condargmat[1].size() == 1 && trueargmat.size() == 2 && trueargmat[1].size() == 1 && falseargmat.size() == 2 && falseargmat[1].size() == 1)
    {
        double* condval = condargmat[1][0].getvalues();
        double* trueval = trueargmat[1][0].getvalues();
        double* falseval = falseargmat[1][0].getvalues();

        for (int i = 0; i < condargmat[1][0].count(); i++)
        {
            if (condval[i] < 0)
                trueval[i] = falseval[i]; 
        }
        
        if (reuse && universe::isreuseallowed)
            universe::setprecomputed(shared_from_this(), trueargmat);
        
        return trueargmat;
    }

    logs log;
    log.msg() << "Error in 'opcondition' object: without FFT the conditional operation can only be computed for constant (harmonic 1) operations" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat opcondition::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat condargmat = mycond->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    densemat trueargmat = mytrue->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    densemat falseargmat = myfalse->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);

    double* condval = condargmat.getvalues();
    double* trueval = trueargmat.getvalues();
    double* falseval = falseargmat.getvalues();

    for (int i = 0; i < condargmat.count(); i++)
    {
        if (condval[i] < 0)
            trueval[i] = falseval[i]; 
    }
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), trueargmat);
    
    return trueargmat;
}

std::shared_ptr<operation> opcondition::simplify(std::vector<int> disjregs)
{
    mycond = mycond->simplify(disjregs);
    mytrue = mytrue->simplify(disjregs);
    myfalse = myfalse->simplify(disjregs);
    
    if (mycond->isconstant() && mycond->getvalue() >= 0)
        return mytrue;
    if (mycond->isconstant() && mycond->getvalue() < 0)
        return myfalse;
    return shared_from_this();
}

std::shared_ptr<operation> opcondition::copy(void)
{
    std::shared_ptr<opcondition> op(new opcondition(mycond,mytrue,myfalse));
    *op = *this;
    op->reuse = false;
    return op;
}

double opcondition::evaluate(void)
{
    if (mycond->evaluate() < 0)
        return myfalse->evaluate();
    else
        return mytrue->evaluate();
}

std::vector<double> opcondition::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaldcond = mycond->evaluate(xcoords, ycoords, zcoords);
    std::vector<double> evaldtrue = mytrue->evaluate(xcoords, ycoords, zcoords);
    std::vector<double> evaldfalse = myfalse->evaluate(xcoords, ycoords, zcoords);

    for (int i = 0; i < evaldcond.size(); i++)
    {
        if (evaldcond[i] < 0)
            evaldtrue[i] = evaldfalse[i];
    }
    return evaldtrue;
}

void opcondition::print(void)
{
    std::cout << "(";
    mycond->print();
    std::cout << " ? ";
    mytrue->print();
    std::cout << ", ";
    myfalse->print();
    std::cout << ")";
}
