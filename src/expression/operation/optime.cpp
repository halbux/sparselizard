#include "optime.h"


std::vector<std::vector<densematrix>> optime::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    densematrix output(elemselect.countinselection(), evaluationcoordinates.size()/3, universe::currenttimestep);
    
    // This can only be on the cos0 harmonic:
    return {{},{output}};
}

densematrix optime::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    std::cout << "Error in 'optime' object: time variable 't' is not supported (non-periodic)" << std::endl;
    abort();
}

std::shared_ptr<operation> optime::copy(void)
{
    std::shared_ptr<optime> op(new optime);
    *op = *this;
    op->reuse = false;
    return op;
}

void optime::print(void)
{
    std::cout << "t";
}
