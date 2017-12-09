#include "oprefelemmeasure.h"


std::vector<std::vector<densematrix>> oprefelemmeasure::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
	int elemtypenum = elemselect.getelementtypenumber();
	
	// Get the measure of the reference element:
	element myelement(elemselect.getelementtypenumber());
	double refelemmeasure = myelement.measurereferenceelement();
	
	densematrix output(elemselect.countinselection(), evaluationcoordinates.size()/3, refelemmeasure);
    
    // On the cos0 harmonic:
	return {{},{output}};
}

densematrix oprefelemmeasure::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
	int elemtypenum = elemselect.getelementtypenumber();
	
	// Get the measure of the reference element:
	element myelement(elemselect.getelementtypenumber());
	double refelemmeasure = myelement.measurereferenceelement();
	
    densematrix output(numtimeevals, elemselect.countinselection() * evaluationcoordinates.size()/3, refelemmeasure);
    
    return output;
}

std::shared_ptr<operation> oprefelemmeasure::copy(void)
{
    std::shared_ptr<oprefelemmeasure> op(new oprefelemmeasure);
    *op = *this;
    return op;
}

void oprefelemmeasure::print(void)
{
    std::cout << "refelemmeasure()";
}
