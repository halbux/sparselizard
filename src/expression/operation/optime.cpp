#include "optime.h"


std::vector<std::vector<densematrix>> optime::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    if (universe::fundamentalfrequency > 0)
    {
        std::cout << "Error in 'optime' object: the time variable 't' cannot be computed without FFT in harmonic domain" << std::endl;
        abort();
    }

    densematrix output(elemselect.countinselection(), evaluationcoordinates.size()/3, universe::currenttimestep);
    return {{},{output}};
}

densematrix optime::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }
    
    int ncols = elemselect.countinselection() * evaluationcoordinates.size()/3;
    densematrix output(numtimeevals, ncols);
    double* outptr = output.getvalues();
    
    double period = 1.0/universe::getfundamentalfrequency();
    double dt = period/numtimeevals;
    
    for (int i = 0; i < numtimeevals; i++)
    {
        double tval = dt*i;
        for (int j = 0; j < ncols; j++)
            outptr[i*ncols+j] = tval;
    }
            
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), output);
    
    return output;
}

std::shared_ptr<operation> optime::copy(void)
{
    std::shared_ptr<optime> op(new optime);
    *op = *this;
    op->reuse = false;
    return op;
}

double optime::evaluate(void)
{
    if (universe::fundamentalfrequency <= 0)
        return universe::currenttimestep;
    else
    {
        std::cout << "Error in 'optime' object: the time variable 't' cannot be evaluated in harmonic domain" << std::endl;
        abort();
    }
}

std::vector<double> optime::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    if (universe::fundamentalfrequency <= 0)
        return std::vector<double>(xcoords.size(), universe::currenttimestep);
    else
    {
        std::cout << "Error in 'optime' object: the time variable 't' cannot be evaluated in harmonic domain" << std::endl;
        abort();
    }
}

void optime::print(void)
{
    std::cout << "t";
}
