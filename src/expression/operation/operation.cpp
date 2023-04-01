#include "operation.h"


std::vector<std::vector<densemat>> operation::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot interpolate the operation" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

densemat operation::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot interpolate the multiharmonic operation" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::vector<densemat>> operation::interpolate(int kietaphiderivative, elementselector& elemselect, std::vector<double>& evaluationcoordinates)
{
    logs log;
    log.msg() << "Error in 'operation' object: expression provided for mesh deformation is invalid" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
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

bool operation::isparameterincluded(std::vector<int> disjregs, rawparameter* rp)
{
    std::vector<std::shared_ptr<operation>> arguments = getarguments();

    for (int i = 0; i < arguments.size(); i++)
    {
        if (arguments[i]->isparameterincluded(disjregs, rp))
            return true;
    }
    return false;
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
    logs log;
    log.msg() << "Error in 'operation' object: cannot get the rawparameter pointer" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::shared_ptr<rawport> operation::getportpointer(void)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot get the rawport pointer" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::shared_ptr<rawfield> operation::getfieldpointer(void)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot get the field pointer" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

bool operation::isreused(void)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot call 'isreused' on the operation" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}
        
void operation::setspacederivative(int whichderivative)
{
    logs log;
    log.msg() << "Error in 'operation' object: either you are trying to apply a space derivative to something else than fields, dof() and tf() or the field does not allow this kind of space derivative" << std::endl;
    log.error();
}   

void operation::setkietaphiderivative(int whichderivative)
{
    logs log;
    log.msg() << "Error in 'operation' object: can only apply reference-element space derivatives to fields, dof() and tf()" << std::endl;
    log.error();
}   

void operation::increasetimederivativeorder(int derivativeorder)
{
    logs log;
    log.msg() << "Error in 'operation' object: can only apply time derivatives to ports, fields, dof() and tf()" << std::endl;
    log.error();
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
    logs log;
    log.msg() << "Error in 'operation' object: cannot copy the operation" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

double operation::evaluate(void)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot evaluate the operation to a double" << std::endl;
    log.msg() << "Did you try to evaluate a space-dependent operation?" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<double> operation::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    logs log;
    log.msg() << "Error in 'operation' object: cannot evaluate the operation to a double" << std::endl;
    log.msg() << "Did you try to evaluate a space-dependent operation?" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

double operation::evaluateattime(double tm)
{
    double tmbkp = universe::currenttimestep;
    
    universe::currenttimestep = tm;
    double evaled = evaluate();
    universe::currenttimestep = tmbkp;
    
    return evaled;    
}

