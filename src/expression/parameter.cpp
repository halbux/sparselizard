#include "parameter.h"


parameter::parameter(void)
{ 
    if (universe::myrawmesh == NULL)
    {
        logs log;
        log.msg() << "Error in 'parameter' object: cannot define a parameter before the mesh is loaded" << std::endl;
        log.error();
    }

    rawparamptr = std::shared_ptr<rawparameter>(new rawparameter(1,1));
}

parameter::parameter(int numrows, int numcols)
{ 
    if (universe::myrawmesh == NULL)
    {
        logs log;
        log.msg() << "Error in 'parameter' object: cannot define a parameter before the mesh is loaded" << std::endl;
        log.error();
    }
    
    rawparamptr = std::shared_ptr<rawparameter>(new rawparameter(numrows,numcols));
}

int parameter::countrows(void)
{
    return rawparamptr->countrows();
}

int parameter::countcolumns(void)
{
    return rawparamptr->countcolumns();
}

parameterselectedregion parameter::operator|(int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    
    return parameterselectedregion(rawparamptr, physreg);
}

void parameter::setvalue(int physreg, expression input)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    
    rawparamptr->set(physreg, input);
}

void parameter::print(void)
{
    rawparamptr->print();
}


vec parameter::atbarycenter(int physreg, field onefield)
{ return ((expression)*this).atbarycenter(physreg, onefield); }

std::vector<double> parameter::max(int physreg, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).max(physreg, refinement, xyzrange); }
std::vector<double> parameter::max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).max(physreg, meshdeform, refinement, xyzrange); }
std::vector<double> parameter::min(int physreg, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).min(physreg, refinement, xyzrange); }
std::vector<double> parameter::min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).min(physreg, meshdeform, refinement, xyzrange); }

void parameter::interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{ ((expression)*this).interpolate(physreg, xyzcoord, interpolated, isfound); }
void parameter::interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{ ((expression)*this).interpolate(physreg, meshdeform, xyzcoord, interpolated, isfound); }

std::vector<double> parameter::interpolate(int physreg, const std::vector<double> xyzcoord)
{ return ((expression)*this).interpolate(physreg, xyzcoord); }
std::vector<double> parameter::interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord)
{ return ((expression)*this).interpolate(physreg, meshdeform, xyzcoord); }
        
double parameter::integrate(int physreg, expression meshdeform, int integrationorder) { return ((expression)*this).integrate(physreg, meshdeform, integrationorder); }
double parameter::integrate(int physreg, int integrationorder) { return ((expression)*this).integrate(physreg, integrationorder); }

void parameter::write(int physreg, int numfftharms, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, filename, lagrangeorder); }
void parameter::write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, meshdeform, filename, lagrangeorder); }

void parameter::write(int physreg, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, filename, lagrangeorder, numtimesteps); }
void parameter::write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, meshdeform, filename, lagrangeorder, numtimesteps); }


expression parameter::operator+(void) { return (expression)*this; }
expression parameter::operator-(void) { return -(expression)*this; }

expression parameter::operator+(parameter inputparameter) { return (expression)*this + inputparameter; }
expression parameter::operator-(parameter inputparameter) { return (expression)*this - inputparameter; }
expression parameter::operator*(parameter inputparameter) { return (expression)*this * inputparameter; }
expression parameter::operator/(parameter inputparameter) { return (expression)*this / inputparameter; }

expression parameter::operator+(double val) { return (expression)*this + val; }
expression parameter::operator-(double val) { return (expression)*this - val; }
expression parameter::operator*(double val) { return (expression)*this * val; }
expression parameter::operator/(double val) { return (expression)*this / val; }


expression operator+(double val, parameter inputparameter) { return (expression)val + inputparameter; }
expression operator-(double val, parameter inputparameter) { return (expression)val - inputparameter; }
expression operator*(double val, parameter inputparameter) { return (expression)val * inputparameter; }
expression operator/(double val, parameter inputparameter) { return (expression)val / inputparameter; }

