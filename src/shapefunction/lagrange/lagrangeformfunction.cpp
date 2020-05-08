#include "lagrangeformfunction.h"

void lagrangeformfunction::preparepoly(void)
{
    if (myformfunctionpolynomials.size() > 0)
        return;

    switch (myelementtypenumber)
    {
        case 0:
            myformfunctionpolynomials = lagrangepoint::getformfunctionpolynomials(myorder);
            break;
        case 1:
            myformfunctionpolynomials = lagrangeline::getformfunctionpolynomials(myorder);
            break;
        case 2:
            myformfunctionpolynomials = lagrangetriangle::getformfunctionpolynomials(myorder);
            break;
        case 3:
            myformfunctionpolynomials = lagrangequadrangle::getformfunctionpolynomials(myorder);
            break;
        case 4:
            myformfunctionpolynomials = lagrangetetrahedron::getformfunctionpolynomials(myorder);
            break;
        case 5:
            myformfunctionpolynomials = lagrangehexahedron::getformfunctionpolynomials(myorder);
            break;
        case 6:
            myformfunctionpolynomials = lagrangeprism::getformfunctionpolynomials(myorder);
            break;
        case 7:
            myformfunctionpolynomials = lagrangepyramid::getformfunctionpolynomials(myorder);
            break;
    }
}

void lagrangeformfunction::preparecoords(void)
{
    if (mynodecoordinates.size() > 0)
        return;
        
    switch (myelementtypenumber)
    {
        case 0:
            mynodecoordinates = lagrangepoint::getnodecoordinates(myorder);
            break;
        case 1:
            mynodecoordinates = lagrangeline::getnodecoordinates(myorder);
            break;
        case 2:
            mynodecoordinates = lagrangetriangle::getnodecoordinates(myorder);
            break;
        case 3:
            mynodecoordinates = lagrangequadrangle::getnodecoordinates(myorder);
            break;
        case 4:
            mynodecoordinates = lagrangetetrahedron::getnodecoordinates(myorder);
            break;
        case 5:
            mynodecoordinates = lagrangehexahedron::getnodecoordinates(myorder);
            break;
        case 6:
            mynodecoordinates = lagrangeprism::getnodecoordinates(myorder);
            break;
        case 7:
            mynodecoordinates = lagrangepyramid::getnodecoordinates(myorder);
            break;
    }
    
    // To avoid problems at the y axis in axisymmetric simulations the nodes are brought slightly closer to the barycenter:
    if (universe::isaxisymmetric)
    {
        int numnodes = mynodecoordinates.size()/3;
    
        std::vector<double> curbarycenter(3,0.0);
        for (int i = 0; i < numnodes; i++)
        {
            curbarycenter[0] += mynodecoordinates[3*i+0]/numnodes;
            curbarycenter[1] += mynodecoordinates[3*i+1]/numnodes;
            curbarycenter[2] += mynodecoordinates[3*i+2]/numnodes;
        }
        
        double lambda = 1e-15;
        for (int i = 0; i < numnodes; i++)
        {
            mynodecoordinates[3*i+0] = (1.0-lambda)*mynodecoordinates[3*i+0] + lambda*curbarycenter[0];
            mynodecoordinates[3*i+1] = (1.0-lambda)*mynodecoordinates[3*i+1] + lambda*curbarycenter[1];
            mynodecoordinates[3*i+2] = (1.0-lambda)*mynodecoordinates[3*i+2] + lambda*curbarycenter[2];
        }
    }
}

lagrangeformfunction::lagrangeformfunction(int elementtypenumber, int order, const std::vector<double> evaluationpoints)
{    
    myorder = order;
    myelementtypenumber = elementtypenumber;
    myevaluationpoints = evaluationpoints;
}

densematrix lagrangeformfunction::getderivative(int whichderivative)
{
    preparepoly();
    
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

std::vector<double> lagrangeformfunction::getnodecoordinates(void)
{
    preparecoords();
    return mynodecoordinates;
}

std::vector<polynomial> lagrangeformfunction::getformfunctionpolynomials(void)
{
    preparepoly();
    return myformfunctionpolynomials;
}

void lagrangeformfunction::print(void)
{
    preparepoly();
    preparecoords();
    
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

