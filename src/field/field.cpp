#include "field.h"


field::field(std::string fieldtypename) { rawfieldptr = shared_ptr<rawfield>(new rawfield(fieldtypename, {1}, false)); }
field::field(std::string fieldtypename, const std::vector<int> harmonicnumbers) 
{ 
    // Make sure all harmonic numbers are positive and non zero:
    for (int i = 0; i < harmonicnumbers.size(); i++)
    {
        if (harmonicnumbers[i] <= 0)
        {
            std::cout << "Error in 'field' object: cannot use negative or zero harmonic number " << harmonicnumbers[i] << std::endl;
            abort();
        }
    }
    if (harmonicnumbers.size() > 0)
        rawfieldptr = shared_ptr<rawfield>(new rawfield(fieldtypename, harmonicnumbers, true)); 
    else
    {
        std::cout << "Error in 'field' object: provided an empty harmonic number list" << std::endl;
        abort();
    } 
}

int field::countcomponents(void) { return rawfieldptr->countcomponents(); }

std::vector<int> field::getharmonics(void) { return rawfieldptr->getharmonics(); }
void field::printharmonics(void) { rawfieldptr->printharmonics(); }

void field::setname(std::string name) { rawfieldptr->setname(name); }
void field::print(void) { rawfieldptr->print(); }

void field::setorder(int physreg, int interpolorder) 
{ 
    if (interpolorder <= 0)
    {
        std::cout << "Error in 'field' object: cannot use negative or zero interpolation order " << interpolorder << std::endl;
        abort();   
    }
    rawfieldptr->setorder(physreg, interpolorder); 
}

void field::setvalue(int physreg, expression input, int extraintegrationdegree) { rawfieldptr->setvalue(physreg, input, extraintegrationdegree); }
void field::setvalue(int physreg, expression input, std::string option) { rawfieldptr->setvalue(physreg, input, option); };
void field::setvalue(int physreg) { rawfieldptr->setvalue(physreg); }

void field::setconstraint(int physreg, expression input, int extraintegrationdegree) { rawfieldptr->setconstraint(physreg, input, extraintegrationdegree); }
void field::setconstraint(int physreg) { rawfieldptr->setconstraint(physreg); }

void field::setdata(int physreg, vectorfieldselect myvec) { rawfieldptr->setdata(physreg, myvec); }
void field::setdata(int physreg, vec myvec) 
{  
    field thisfield = *this;
    setdata(physreg, myvec|thisfield); 
}

field field::comp(int component) 
{ 
    if (component < 0 || component > 2)
    {
        std::cout << "Error in 'field' object: cannot use component number " << component << " (only 0, 1 and 2 are allowed)" << std::endl;
        abort();
    }
    return field(rawfieldptr->comp(component)); 
}

field field::harmonic(const std::vector<int> harmonicnumbers)
{
    if (harmonicnumbers.size() == 0)
    {
        std::cout << "Error in 'field' object: no harmonics provided to the .harmonic function" << std::endl;
        abort();
    }    
    // Make sure all harmonic numbers are positive and non zero:
    for (int i = 0; i < harmonicnumbers.size(); i++)
    {
        if (harmonicnumbers[i] <= 0)
        {
            std::cout << "Error in 'field' object: cannot use negative or zero harmonic number " << harmonicnumbers[i] << std::endl;
            abort();
        }
    }
    return field(rawfieldptr->harmonic(harmonicnumbers));
}




double field::integrate(int physreg, expression meshdeform, int integrationorder) { return ((expression)*this).integrate(physreg, meshdeform, integrationorder); }
double field::integrate(int physreg, int integrationorder) { return ((expression)*this).integrate(physreg, integrationorder); }

void field::write(int physreg, int numfftharms, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, filename, lagrangeorder); }
void field::write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, meshdeform, filename, lagrangeorder); }

void field::write(int physreg, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, filename, lagrangeorder, numtimesteps); }
void field::write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, meshdeform, filename, lagrangeorder, numtimesteps); }


expression field::operator+(void) { return (expression)*this; }
expression field::operator-(void) { return -(expression)*this; }

expression field::operator+(field inputfield) { return (expression)*this + inputfield; }
expression field::operator-(field inputfield) { return (expression)*this - inputfield; }
expression field::operator*(field inputfield) { return (expression)*this * inputfield; }
expression field::operator/(field inputfield) { return (expression)*this / inputfield; }

expression field::operator+(double val) { return (expression)*this + val; }
expression field::operator-(double val) { return (expression)*this - val; }
expression field::operator*(double val) { return (expression)*this * val; }
expression field::operator/(double val) { return (expression)*this / val; }

expression field::operator+(parameter& param) { return (expression)*this + param; }
expression field::operator-(parameter& param) { return (expression)*this - param; }
expression field::operator*(parameter& param) { return (expression)*this * param; }
expression field::operator/(parameter& param) { return (expression)*this / param; }       
 

expression operator+(double val, field inputfield) { return inputfield+val; }
expression operator-(double val, field inputfield) { return -inputfield+val; }
expression operator*(double val, field inputfield) { return inputfield*val; }
expression operator/(double val, field inputfield) { return ( (expression)val )/( (expression)inputfield ); }

expression operator+(parameter& param, field inputfield) { return inputfield+param; }
expression operator-(parameter& param, field inputfield) { return -inputfield+param; }
expression operator*(parameter& param, field inputfield) { return inputfield*param; }
expression operator/(parameter& param, field inputfield) { return ( (expression)param ) / ( (expression)inputfield ); }


