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
field::field(std::string fieldtypename, spanningtree spantree)
{
	rawfieldptr = shared_ptr<rawfield>(new rawfield(fieldtypename, {1}, false));
	rawfieldptr->setspanningtree(spantree);
}

field::field(std::string fieldtypename, const std::vector<int> harmonicnumbers, spanningtree spantree)
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
	rawfieldptr->setspanningtree(spantree);
}

int field::countcomponents(void) { return rawfieldptr->countcomponents(); }

std::vector<int> field::getharmonics(void) { return rawfieldptr->getharmonics(); }
void field::printharmonics(void) { rawfieldptr->printharmonics(); }

void field::setname(std::string name) { rawfieldptr->setname(name); }
void field::print(void) { rawfieldptr->print(); }

void field::setorder(int physreg, int interpolorder) 
{ 
    if (interpolorder < 0)
    {
        std::cout << "Error in 'field' object: cannot use negative interpolation order " << interpolorder << std::endl;
        abort();   
    }
    if (interpolorder == 0 && rawfieldptr->gettypename() != "hcurl")
    {
        std::cout << "Error in 'field' object: cannot use interpolation order 0 for shape function " << rawfieldptr->gettypename() << std::endl;
        abort();   
    }
    rawfieldptr->setorder(physreg, interpolorder); 
}

void field::setvalue(int physreg, expression input, int extraintegrationdegree) { rawfieldptr->setvalue(physreg, -1, NULL, input, extraintegrationdegree); }
void field::setvalue(int physreg, expression meshdeform, expression input, int extraintegrationdegree) { rawfieldptr->setvalue(physreg, -1, &meshdeform, input, extraintegrationdegree); }
void field::setvalue(int physreg, int numfftharms, expression input, int extraintegrationdegree) { rawfieldptr->setvalue(physreg, numfftharms, NULL, input, extraintegrationdegree); }
void field::setvalue(int physreg, int numfftharms, expression meshdeform, expression input, int extraintegrationdegree) { rawfieldptr->setvalue(physreg, numfftharms, &meshdeform, input, extraintegrationdegree); }
void field::setvalue(int physreg) { rawfieldptr->setvalue(physreg); }

void field::setconstraint(int physreg, expression input, int extraintegrationdegree) { rawfieldptr->setconstraint(physreg, input, extraintegrationdegree); }
void field::setconstraint(int physreg) { rawfieldptr->setconstraint(physreg); }

void field::setconditionalconstraint(int physreg, expression condexpr, expression valexpr) { rawfieldptr->setconditionalconstraint(physreg, condexpr, valexpr); }

void field::setgauge(int physreg) 
{ 
    if (rawfieldptr->gettypename() != "hcurl")
    {
        std::cout << "Error in 'field' object: cannot gauge shape function " << rawfieldptr->gettypename() << " (only hcurl)" << std::endl;
        abort();   
    }

	rawfieldptr->setgauge(physreg);
}

void field::setdata(int physreg, vectorfieldselect myvec, std::string op) 
{ 
	if (op != "set" && op != "add")
	{
		std::cout << "Error in 'field' object: operation " << op << " is unknown in .setdata (use 'set' or 'add')" << std::endl;
		abort();
	}

	rawfieldptr->setdata(physreg, myvec, op); 
}

void field::setdata(int physreg, vec myvec, std::string op)
{  
    field thisfield = *this;
    setdata(physreg, myvec|thisfield, op); 
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




vec field::atbarycenter(int physreg, field onefield) { return ((expression)*this).atbarycenter(physreg, onefield); }

std::vector<double> field::max(int physreg, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).max(physreg, refinement, xyzrange); }
std::vector<double> field::max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).max(physreg, meshdeform, refinement, xyzrange); }
std::vector<double> field::min(int physreg, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).min(physreg, refinement, xyzrange); }
std::vector<double> field::min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).min(physreg, meshdeform, refinement, xyzrange); }

void field::interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{ ((expression)*this).interpolate(physreg, xyzcoord, interpolated, isfound); }
void field::interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{ ((expression)*this).interpolate(physreg, meshdeform, xyzcoord, interpolated, isfound); }

std::vector<double> field::interpolate(int physreg, const std::vector<double> xyzcoord)
{ return ((expression)*this).interpolate(physreg, xyzcoord); }
std::vector<double> field::interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord)
{ return ((expression)*this).interpolate(physreg, meshdeform, xyzcoord); }

double field::integrate(int physreg, expression meshdeform, int integrationorder) { return ((expression)*this).integrate(physreg, meshdeform, integrationorder); }
double field::integrate(int physreg, int integrationorder) { return ((expression)*this).integrate(physreg, integrationorder); }

void field::write(int physreg, int numfftharms, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, filename, lagrangeorder); }
void field::write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, meshdeform, filename, lagrangeorder); }

void field::write(int physreg, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, filename, lagrangeorder, numtimesteps); }
void field::write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, meshdeform, filename, lagrangeorder, numtimesteps); }

void field::writeraw(int physreg, std::string filename, bool isbinary)
{
    if (rawfieldptr == NULL)
    {
        std::cout << "Error in 'field' object: cannot write raw field data (field undefined)" << std::endl;
        abort();
    }

    if (isbinary == true && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz" || filename.size() >= 8 && filename.substr(filename.size()-7,7) == ".slz.gz"))
        rawfieldptr->writeraw(physreg, filename, isbinary);
    
    if (isbinary == false && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz"))
        rawfieldptr->writeraw(physreg, filename, isbinary);
    
    std::cout << "Error in 'field' object: expected .slz file extension (or .slz.gz for binary format) to write raw field data" << std::endl;
    abort();
}

void field::loadraw(std::string filename, bool isbinary)
{
    if (rawfieldptr == NULL)
    {
        std::cout << "Error in 'field' object: cannot load raw field data (field undefined)" << std::endl;
        abort();
    }

    if (isbinary == true && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz" || filename.size() >= 8 && filename.substr(filename.size()-7,7) == ".slz.gz"))
        rawfieldptr->loadraw(filename, isbinary);
    
    if (isbinary == false && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz"))
        rawfieldptr->loadraw(filename, isbinary);
    
    std::cout << "Error in 'field' object: expected .slz file extension (or .slz.gz for binary format) to load raw field data" << std::endl;
    abort();
}


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


