#include "opcustom.h"


opcustom::opcustom(int outindex, std::vector<densemat> fct(std::vector<densemat>), std::vector<std::shared_ptr<operation>> args)
{
    myoutindex = outindex;
    myargs = args;
    myfunction = fct;
}

opcustom::opcustom(int outindex, std::vector<densemat> fct(std::vector<densemat>, std::vector<field>, elementselector&, std::vector<double>&, expression*), std::vector<std::shared_ptr<operation>> args, std::vector<field> infields)
{
    myoutindex = outindex;
    myargs = args;
    myadvancedfunction = fct;
    myfields = infields;
}
        
std::vector<std::vector<densemat>> opcustom::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available:
    if (universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    int numels = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    std::vector<densemat> fctargs(myargs.size());
    for (int i = 0; i < myargs.size(); i++)
    {
        std::vector<std::vector<densemat>> argmat = myargs[i]->interpolate(elemselect, evaluationcoordinates, meshdeform);
        if (argmat.size() != 2 || argmat[1].size() != 1)
        {
            logs log;
            log.msg() << "Error in 'opcustom' object: without FFT the custom operation can only be computed for constant (harmonic 1) operations" << std::endl;
            log.error();
        }
        fctargs[i] = argmat[1][0];
    }
    
    std::vector<densemat> output;
    if (myfunction != NULL)
    {
        bool wasreuseallowed = universe::isreuseallowed;
        auto storage = universe::backup();
        // Safe call to custom function:
        universe::forbidreuse();
        output = myfunction(fctargs);
        universe::forbidreuse();
        
        universe::restore(storage);
        if (wasreuseallowed)
            universe::allowreuse();
    }
    else
        output = myadvancedfunction(fctargs, myfields, elemselect, evaluationcoordinates, meshdeform);
    
    // Make sure the user provided function returns something valid:
    if (output.size() != myfamily.size())
    {
        logs log;
        log.msg() << "Error in 'opcustom' object: custom function returned " << output.size() << " densemat objects (expected " << myfamily.size() <<  ")" << std::endl;
        log.error();
    }
    for (int i = 0; i < output.size(); i++)
    {
        if (output[i].isdefined() == false)
        {
            logs log;
            log.msg() << "Error in 'opcustom' object: custom function returned an undefined densemat object" << std::endl;
            log.error();
        }
        if (output[i].countrows() != numels || output[i].countcolumns() != numevalpts)
        {
            logs log;
            log.msg() << "Error in 'opcustom' object: custom function returned a " << output[i].countrows() << "x" << output[i].countcolumns() << " densemat (expected " << numels << "x" << numevalpts << ")" << std::endl;
            log.error();
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

densemat opcustom::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available:
    if (universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    int numels = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    std::vector<densemat> fctargs(myargs.size());
    for (int i = 0; i < myargs.size(); i++)
        fctargs[i] = myargs[i]->multiharmonicinterpolate(numtimeevals, elemselect, evaluationcoordinates, meshdeform);
    
    std::vector<densemat> output;
    if (myfunction != NULL)
    {
        bool wasreuseallowed = universe::isreuseallowed;
        auto storage = universe::backup();
        // Safe call to custom function:
        universe::forbidreuse();
        output = myfunction(fctargs);
        universe::forbidreuse();
        
        universe::restore(storage);
        if (wasreuseallowed)
            universe::allowreuse();
    }
    else
        output = myadvancedfunction(fctargs, myfields, elemselect, evaluationcoordinates, meshdeform);
    
    // Make sure the user provided function returns something valid:
    if (output.size() != myfamily.size())
    {
        logs log;
        log.msg() << "Error in 'opcustom' object: custom function returned " << output.size() << " densemat objects (expected " << myfamily.size() <<  ")" << std::endl;
        log.error();
    }
    for (int i = 0; i < output.size(); i++)
    {
        if (output[i].isdefined() == false)
        {
            logs log;
            log.msg() << "Error in 'opcustom' object: custom function returned an undefined densemat object" << std::endl;
            log.error();
        }
        if (output[i].countrows() != numtimeevals || output[i].countcolumns() != numels*numevalpts)
        {
            logs log;
            log.msg() << "Error in 'opcustom' object: custom function returned a " << output[i].countrows() << "x" << output[i].countcolumns() << " densemat (expected " << numtimeevals << "x" << numels*numevalpts << ")" << std::endl;
            log.error();
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
    std::shared_ptr<opcustom> op;
    if (myfunction != NULL)
        op = std::shared_ptr<opcustom>(new opcustom(myoutindex, myfunction, myargs));
    else
        op = std::shared_ptr<opcustom>(new opcustom(myoutindex, myadvancedfunction, myargs, myfields));
    *op = *this;
    return op;
}

void opcustom::print(void)
{
    std::cout << "custom";
}
