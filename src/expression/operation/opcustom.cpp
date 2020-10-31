#include "opcustom.h"


opcustom::opcustom(int outindex, std::vector<densematrix> fct(std::vector<densematrix>), std::vector<std::shared_ptr<operation>> args)
{
    myoutindex = outindex;
    myargs = args;
    myfunction = fct;
}
        
std::vector<std::vector<densematrix>> opcustom::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available:
    if (universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    int numels = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    std::vector<densematrix> fctargs(myargs.size());
    for (int i = 0; i < myargs.size(); i++)
    {
        std::vector<std::vector<densematrix>> argmat = myargs[i]->interpolate(elemselect, evaluationcoordinates, meshdeform);
        if (argmat.size() != 2 || argmat[1].size() != 1)
        {
            std::cout << "Error in 'opcustom' object: without FFT the custom operation can only be computed for constant (harmonic 1) operations" << std::endl;
            abort();
        }
        fctargs[i] = argmat[1][0];
    }
    
    std::vector<densematrix> output = myfunction(fctargs);
    
    // Make sure the user provided function returns something valid:
    if (output.size() != myfamily.size())
    {
        std::cout << "Error in 'opcustom' object: custom function returned " << output.size() << " densematrix objects (expected " << myfamily.size() <<  ")" << std::endl;
        abort();
    }
    for (int i = 0; i < output.size(); i++)
    {
        if (output[i].isdefined() == false)
        {
            std::cout << "Error in 'opcustom' object: custom function returned an undefined densematrix object" << std::endl;
            abort();
        }
        if (output[i].countrows() != numels || output[i].countcolumns() != numevalpts)
        {
            std::cout << "Error in 'opcustom' object: custom function returned a " << output[i].countrows() << "x" << output[i].countcolumns() << " densematrix (expected " << numels << "x" << numevalpts << ")" << std::endl;
            abort();
        }
    }
    
    if (universe::isreuseallowed)
    {
        for (int i = 0; i < myfamily.size(); i++)
        {
            if (myfamily[i].expired() == false)
                universe::setprecomputed(myfamily[i].lock(), {{},{output[i]}});
        }
    }
    
    return {{},{output[myoutindex]}};
}

densematrix opcustom::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available:
    if (universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    int numels = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    std::vector<densematrix> fctargs(myargs.size());
    for (int i = 0; i < myargs.size(); i++)
        fctargs[i] = myargs[i]->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    
    std::vector<densematrix> output = myfunction(fctargs);
    
    // Make sure the user provided function returns something valid:
    if (output.size() != myfamily.size())
    {
        std::cout << "Error in 'opcustom' object: custom function returned " << output.size() << " densematrix objects (expected " << myfamily.size() <<  ")" << std::endl;
        abort();
    }
    for (int i = 0; i < output.size(); i++)
    {
        if (output[i].isdefined() == false)
        {
            std::cout << "Error in 'opcustom' object: custom function returned an undefined densematrix object" << std::endl;
            abort();
        }
        if (output[i].countrows() != numtimeevals || output[i].countcolumns() != numels*numevalpts)
        {
            std::cout << "Error in 'opcustom' object: custom function returned a " << output[i].countrows() << "x" << output[i].countcolumns() << " densematrix (expected " << numtimeevals << "x" << numels*numevalpts << ")" << std::endl;
            abort();
        }
    }
    
    if (universe::isreuseallowed)
    {
        for (int i = 0; i < myfamily.size(); i++)
        {
            if (myfamily[i].expired() == false)
                universe::setprecomputedfft(myfamily[i].lock(), output[i]);
        }
    }
    
    return output[myoutindex];
}

std::shared_ptr<operation> opcustom::simplify(std::vector<int> disjregs)
{
    for (int i = 0; i < myargs.size(); i++)
        myargs[i] = myargs[i]->simplify(disjregs);
    
    return shared_from_this();
}

std::shared_ptr<operation> opcustom::copy(void)
{
    std::shared_ptr<opcustom> op(new opcustom(myoutindex, myfunction, myargs));
    *op = *this;
    return op;
}

void opcustom::print(void)
{
    std::cout << "custom";
}
