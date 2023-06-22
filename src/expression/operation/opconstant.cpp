#include "opconstant.h"


std::vector<std::vector<densemat>> opconstant::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    densemat output(elemselect.countinselection(), evaluationcoordinates.size()/3, constantvalue);
    
    // The constant is on the cos0 harmonic:
    return {{},{output}};
}

densemat opconstant::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    densemat output(numtimeevals, elemselect.countinselection() * evaluationcoordinates.size()/3, constantvalue);
    
    return output;
}

std::shared_ptr<operation> opconstant::copy(void)
{
    std::shared_ptr<opconstant> op(new opconstant(constantvalue));
    *op = *this;
    op->reuse = false;
    return op;
}

std::vector<double> opconstant::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::vector<double> evaluated(xcoords.size(), constantvalue);
    return evaluated;
}

void opconstant::print(void) { std::cout << constantvalue; }
