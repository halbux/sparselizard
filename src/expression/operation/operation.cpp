#include "operation.h"


std::vector<std::vector<densemat>> operation::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    std::cout << "Error in 'operation' object: cannot interpolate the operation" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

densemat operation::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    std::cout << "Error in 'operation' object: cannot interpolate the multiharmonic operation" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

std::vector<std::vector<densemat>> operation::interpolate(int kietaphiderivative, elementselector& elemselect, std::vector<double>& evaluationcoordinates)
{
    std::cout << "Error in 'operation' object: expression provided for mesh deformation is invalid" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

bool operation::iszero(void) { return (isconstant() && getvalue() == 0); }

bool operation::isdofincluded(void)
{ 
    std::vector<std::shared_ptr<operation>> arguments = getarguments();

    for (int i = 0; i < arguments.size(); i++)
    {
        if (arguments[i]->isdofincluded())
            return true;
    }
    return isdof();
}

bool operation::istfincluded(void)
{ 
    std::vector<std::shared_ptr<operation>> arguments = getarguments();
    
    for (int i = 0; i < arguments.size(); i++)
    {
        if (arguments[i]->istfincluded())
            return true;
    }
    return istf();
}

bool operation::isportincluded(void)
{ 
    std::vector<std::shared_ptr<operation>> arguments = getarguments();

    for (int i = 0; i < arguments.size(); i++)
    {
        if (arguments[i]->isportincluded())
            return true;
    }
    return isport();
}

bool operation::isharmonicone(std::vector<int> disjregs)
{
    std::vector<std::shared_ptr<operation>> arguments = getarguments();
    
    for (int i = 0; i < arguments.size(); i++)
    {
        if (arguments[i]->isharmonicone(disjregs) == false)
            return false;
    }
    return true;
}

std::shared_ptr<rawparameter> operation::getparameterpointer(void)
{
    std::cout << "Error in 'operation' object: cannot get the rawparameter pointer" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

std::shared_ptr<rawport> operation::getportpointer(void)
{
    std::cout << "Error in 'operation' object: cannot get the rawport pointer" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

std::shared_ptr<rawfield> operation::getfieldpointer(void)
{
    std::cout << "Error in 'operation' object: cannot get the field pointer" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

bool operation::isreused(void)
{
    std::cout << "Error in 'operation' object: cannot call 'isreused' on the operation" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}
        
void operation::setspacederivative(int whichderivative)
{
    std::cout << "Error in 'operation' object: either you are trying to apply a space derivative to something else than fields, dof() and tf() or the field does not allow this kind of space derivative" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}   

void operation::setkietaphiderivative(int whichderivative)
{
    std::cout << "Error in 'operation' object: can only apply reference-element space derivatives to fields, dof() and tf()" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}   

void operation::increasetimederivativeorder(int derivativeorder)
{
    std::cout << "Error in 'operation' object: can only apply time derivatives to ports, fields, dof() and tf()" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}   

bool operation::isvalueorientationdependent(std::vector<int> disjregs)
{
    std::vector<std::shared_ptr<operation>> arguments = getarguments();
    
    if (arguments.size() > 0)
    {
        for (int i = 0; i < arguments.size(); i++)
        {
            if (arguments[i]->isvalueorientationdependent(disjregs))
                return true;
        }
    }
    return false;
}

std::shared_ptr<operation> operation::copy(void)
{
    std::cout << "Error in 'operation' object: cannot copy the operation" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

double operation::evaluate(void)
{
    std::cout << "Error in 'operation' object: cannot evaluate the operation to a double" << std::endl;
    std::cout << "Did you try to evaluate a space-dependent operation?" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

std::vector<double> operation::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    std::cout << "Error in 'operation' object: cannot evaluate the operation to a double" << std::endl;
    std::cout << "Did you try to evaluate a space-dependent operation?" << std::endl;
    std::cout << "Operation was:" << std::endl;
    this->print();
    std::cout << std::endl;
    abort();
}

double operation::evaluateattime(double tm)
{
    double tmbkp = universe::currenttimestep;
    
    universe::currenttimestep = tm;
    double evaled = evaluate();
    universe::currenttimestep = tmbkp;
    
    return evaled;    
}

