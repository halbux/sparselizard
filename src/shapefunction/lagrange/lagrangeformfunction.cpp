#include "lagrangeformfunction.h"

lagrangeformfunction::lagrangeformfunction(int elementtypenumber, int order, const std::vector<double> evaluationpoints)
{    
    myorder = order;
    myelementtypenumber = elementtypenumber;
	myevaluationpoints = evaluationpoints;
    
    switch (myelementtypenumber)
    {
        case 0:
            mynodecoordinates = lagrangepoint::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangepoint::getformfunctionpolynomials(myorder);
            break;
        case 1:
            mynodecoordinates = lagrangeline::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangeline::getformfunctionpolynomials(myorder);
            break;
        case 2:
            mynodecoordinates = lagrangetriangle::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangetriangle::getformfunctionpolynomials(myorder);
            break;
        case 3:
            mynodecoordinates = lagrangequadrangle::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangequadrangle::getformfunctionpolynomials(myorder);
            break;
        case 4:
            mynodecoordinates = lagrangetetrahedron::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangetetrahedron::getformfunctionpolynomials(myorder);
            break;
        case 5:
            mynodecoordinates = lagrangehexahedron::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangehexahedron::getformfunctionpolynomials(myorder);
            break;
        case 6:
            mynodecoordinates = lagrangeprism::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangeprism::getformfunctionpolynomials(myorder);
            break;
        case 7:
            mynodecoordinates = lagrangepyramid::getnodecoordinates(myorder);
            myformfunctionpolynomials = lagrangepyramid::getformfunctionpolynomials(myorder);
            break;
    }
}

densematrix lagrangeformfunction::getderivative(int whichderivative)
{
    // In case it has already been computed:
    if (evaluated[whichderivative].isdefined())
        return evaluated[whichderivative].copy();
    
    element myelement(myelementtypenumber, myorder);
    int numberofformfunctions = myelement.countcurvednodes();
    int numberofevaluationpoints = myevaluationpoints.size()/3;
    
    evaluated[whichderivative] = densematrix(numberofformfunctions, numberofevaluationpoints);
    
    // Fill the dense matrices with the polynomial value:
    for (int i = 0; i < numberofformfunctions; i++)
        evaluated[whichderivative].setrow(i, myformfunctionpolynomials[i].evalat(myevaluationpoints, whichderivative));

    // Return a copy to make sure it is not changed.
    return evaluated[whichderivative].copy();
}

void lagrangeformfunction::print(void)
{
    element myelement(myelementtypenumber, myorder);
    int numberofformfunctions = myelement.countcurvednodes();
    
    std::cout << "Printing the " << numberofformfunctions << " order " << myorder << " Lagrange form functions for element " << myelement.gettypename() << ":" << std::endl << std::endl;
    int oldprecision = std::cout.precision();
    std::cout.precision(17);
    for (int i = 0; i < numberofformfunctions; i++)
    {
        std::cout << "Number " << i << " on node with reference coordinates (ki eta phi) = (";
        std::cout << std::setw(26) << std::left << mynodecoordinates[3*i+0] << std::setw(26) << std::left << mynodecoordinates[3*i+1] << std::setw(19) << std::left << mynodecoordinates[3*i+2];
        std::cout << "):" << std::endl;
        myformfunctionpolynomials[i].print();
    }
    std::cout.precision(oldprecision);
}

