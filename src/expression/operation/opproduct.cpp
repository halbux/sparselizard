#include "opproduct.h"


std::vector<std::vector<densemat>> opproduct::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    std::vector<std::vector<densemat>> product = productterms[0]->interpolate(elemselect, evaluationcoordinates, meshdeform);
    bool isonlycos0 = (product.size() == 2 && product[1].size() == 1);
    
    // Multiply 'product' by all remaining terms:
    for (int i = 1; i < productterms.size(); i++)
    {
        std::vector<std::vector<densemat>> currentterm = productterms[i]->interpolate(elemselect, evaluationcoordinates, meshdeform);        
        bool iscurrenttermonlycos0 = (currentterm.size() == 2 && currentterm[1].size() == 1);

        // Shortcut for cos0 * cos0:
        if (isonlycos0 && iscurrenttermonlycos0)
        {
            product[1][0].multiplyelementwise(currentterm[1][0]);
            continue;
        }
        if (isonlycos0 && not(iscurrenttermonlycos0))
        {
            for (int harm = 1; harm < currentterm.size(); harm++)
            {
                if (currentterm[harm].size() == 1)
                    currentterm[harm][0].multiplyelementwise(product[1][0]);
            }
            product = currentterm;
            isonlycos0 = false;
            continue;
        }
        if (not(isonlycos0) && iscurrenttermonlycos0)
        {
            for (int harm = 1; harm < product.size(); harm++)
            {
                if (product[harm].size() == 1)
                    product[harm][0].multiplyelementwise(currentterm[1][0]);
            }
            continue;
        }
        
        // General case:
        if (not(isonlycos0) && not(iscurrenttermonlycos0))
        {
            std::vector<std::vector<densemat>> tempproduct = {};
        
            // Loop on all product harmonics:
            for (int pharm = 0; pharm < product.size(); pharm++)
            {
                if (product[pharm].size() == 0)
                    continue;
                
                // Loop on all harmonics of the current term:
                for (int charm = 0; charm < currentterm.size(); charm++)
                {
                    if (currentterm[charm].size() == 0)
                        continue;
                    
                    // Precompute the product of the current set of harmonics:
                    densemat curprod = product[pharm][0].copy();
                    curprod.multiplyelementwise(currentterm[charm][0]);
                
                    std::vector<std::pair<int, double>> harmsofproduct = harmonic::getproduct(pharm, charm);

                    for (int p = 0; p < harmsofproduct.size(); p++)
                    {
                        int currentharm = harmsofproduct[p].first;
                        // currentharmcoef can be + or - 0.5 or 1.
                        double currentharmcoef = harmsofproduct[p].second;
                        
                        if (currentharm == 0)
                            continue;
                        
                        if (tempproduct.size() < currentharm+1)
                            tempproduct.resize(currentharm+1);
                            
                        if (tempproduct[currentharm].size() == 1)
                            tempproduct[currentharm][0].addproduct(currentharmcoef, curprod);
                        else
                            tempproduct[currentharm] = {curprod.getproduct(currentharmcoef)};
                    }
                }
            }
            product = tempproduct;
        }
    }
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), product);
    
    return product;
}

densemat opproduct::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    densemat output = productterms[0]->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    
    for (int i = 1; i < productterms.size(); i++)
        output.multiplyelementwise(productterms[i]->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform));
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> opproduct::expand(void)
{
    // Expand all factors including a dof or tf. It is not required 
    // to expand the other ones for the formulations.
    for (int i = 0; i < productterms.size(); i++)
    {
        // We only want to expand operations that include a dof(), tf() or port.
        if (productterms[i]->isdofincluded() || productterms[i]->istfincluded() || productterms[i]->isportincluded())
            productterms[i] = productterms[i]->expand();
    }
    
    // Move all sub products to this product operation:
    group();
        
    // Everything is first put in form of a sum, possibly with a single sum term.
    // If there is no dof, tf or port this leads to treating a sum as a monolithic block.
    std::vector<std::shared_ptr<operation>> prodtrms(productterms.size());
    for (int i = 0; i < productterms.size(); i++)
    {
        prodtrms[i] = productterms[i];
        if (not(productterms[i]->issum()) || not(productterms[i]->isdofincluded()) && not(productterms[i]->istfincluded()) && not(productterms[i]->isportincluded()))
        {
            std::shared_ptr<opsum> op(new opsum);
            op->addterm(productterms[i]);
            prodtrms[i] = op;
        }
    }
    
    std::shared_ptr<operation> expanded = prodtrms[0];
    
    // Multiply 'expanded' by all other factors, one at a time:
    for (int i = 1; i < prodtrms.size(); i++)
    {
        // Use the sum*sum expansion rule since we have only sums:
        std::shared_ptr<opsum> opersum(new opsum);
        
        for (int j = 0; j < expanded->count(); j++)
        {
            for (int k = 0; k < prodtrms[i]->count(); k++)
            {
                std::shared_ptr<opproduct> prod(new opproduct);

                prod->multiplybyterm(expanded->getargument(j));
                prod->multiplybyterm(prodtrms[i]->getargument(k));
                
                prod->group();
                
                opersum->addterm(prod);
            }
        }
        expanded = opersum;
    }
    // If there is a single sum term:
    if (expanded->count() <= 1)
        expanded = expanded->getargument(0);
        
    return expanded;
}

void opproduct::group(void)
{
    std::vector<std::shared_ptr<operation>> groupedproductterms = {};
    
    for (int i = 0; i < productterms.size(); i++)
    {
        productterms[i]->group();

        // We want to keep the operations to reuse untouched!
        if (productterms[i]->isproduct() && productterms[i]->isreused() == false)
        {
            for (int j = 0; j < productterms[i]->count(); j++)
                groupedproductterms.push_back(productterms[i]->getargument(j));
        }
        else
            groupedproductterms.push_back(productterms[i]);
    }
    productterms = groupedproductterms;
}

std::shared_ptr<operation> opproduct::simplify(std::vector<int> disjregs)
{
    for (int i = 0; i < count(); i++)
        productterms[i] = productterms[i]->simplify(disjregs);
    
    group();

    double multiplied = 1;
    int numproductterms = count();
    for (int i = numproductterms-1; i >= 0; i--)
    {
        if (productterms[i]->isconstant())
        {
            // Any 0 product term leads to a 0 output:
            double currentval = productterms[i]->getvalue();
            if (currentval == 0)
            {
                // We do not want to invalidate this operation and thus add the constant term:
                productterms = {std::shared_ptr<operation>(new opconstant(0))};
                return productterms[0];
            }
            
            multiplied *= currentval;
            removeterm(i);
        }
    }
    // If there were not only constants in the product:
    if (count() > 0)
    {
        if (multiplied != 1)
            multiplybyterm(std::shared_ptr<operation>(new opconstant(multiplied)));
    }
    else
    {
        // We do not want to invalidate this operation and thus add the constant term:
        productterms = {std::shared_ptr<operation>(new opconstant(multiplied))};
        return productterms[0];
    }
    
    return shared_from_this();
}

std::shared_ptr<operation> opproduct::copy(void)
{
    std::shared_ptr<opproduct> op(new opproduct);
    *op = *this;
    op->reuse = false;
    return op;
}

std::vector<double> opproduct::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluated(xcoords.size(), 1);
    for (int i = 0; i < productterms.size(); i++)
    {
        std::vector<double> current = productterms[i]->evaluate(xcoords, ycoords, zcoords);
        for (int j = 0; j < xcoords.size(); j++)
            evaluated[j] *= current[j];
    }
    return evaluated;
}

void opproduct::print(void)
{
    for (int i = 0; i < productterms.size(); i++)
    {
        if (i > 0) 
            std::cout << " * ";
        productterms[i]->print();
    }
}
