#include "field.h"


void field::errorifpointerisnull(void)
{
    if (rawfieldptr == NULL)
    {
        std::cout << "Error in 'field' object: cannot perform the operation (field is undefined)" << std::endl;
        abort();
    }
}

field::field(std::string fieldtypename) { rawfieldptr = std::shared_ptr<rawfield>(new rawfield(fieldtypename, {1}, false)); }
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
        rawfieldptr = std::shared_ptr<rawfield>(new rawfield(fieldtypename, harmonicnumbers, true)); 
    else
    {
        std::cout << "Error in 'field' object: provided an empty harmonic number list" << std::endl;
        abort();
    } 
}
field::field(std::string fieldtypename, spanningtree spantree)
{
    rawfieldptr = std::shared_ptr<rawfield>(new rawfield(fieldtypename, {1}, false));
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
        rawfieldptr = std::shared_ptr<rawfield>(new rawfield(fieldtypename, harmonicnumbers, true)); 
    else
    {
        std::cout << "Error in 'field' object: provided an empty harmonic number list" << std::endl;
        abort();
    } 
    rawfieldptr->setspanningtree(spantree);
}

int field::countcomponents(void) { errorifpointerisnull(); return rawfieldptr->countcomponents(); }

std::vector<int> field::getharmonics(void) { errorifpointerisnull(); return rawfieldptr->getharmonics(); }
void field::printharmonics(void) { errorifpointerisnull(); rawfieldptr->printharmonics(); }

void field::setname(std::string name) { errorifpointerisnull(); rawfieldptr->setname(name); }
void field::print(void) { errorifpointerisnull(); rawfieldptr->print(); }

void field::setorder(int physreg, int interpolorder) 
{ 
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    
    if (interpolorder < 0)
    {
        std::cout << "Error in 'field' object: cannot use negative interpolation order " << interpolorder << std::endl;
        abort();   
    }
    hierarchicalformfunction hff;
    if (interpolorder < hff.getminorder(rawfieldptr->gettypename()))
    {
        std::cout << "Error in 'field' object: cannot use interpolation order " << interpolorder << " for shape function " << rawfieldptr->gettypename() << std::endl;
        abort();   
    }
    rawfieldptr->setorder(physreg, interpolorder); 
}

void field::setorder(expression criterion, int loworder, int highorder)
{
    errorifpointerisnull();
    
    std::string tn = rawfieldptr->gettypename();
    
    if (not(criterion.isscalar()))
    {
        std::cout << "Error in 'field' object: expected a scalar criterion for p-adaptivity" << std::endl;
        abort();   
    }
    // The criterion cannot be multiharmonic:
    std::vector<int> alldisjregs((universe::mymesh->getdisjointregions())->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(criterion.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'field' object: cannot have a multiharmonic criterion for p-adaptivity" << std::endl;
        abort();
    }
    
    if (loworder < 0)
    {
        std::cout << "Error in 'field' object: in 'setorder' cannot use negative minimum interpolation order " << loworder << std::endl;
        abort();   
    }
    hierarchicalformfunction hff;
    if (loworder < hff.getminorder(tn))
    {
        std::cout << "Error in 'field' object: in 'setorder' cannot use interpolation order " << loworder << " for shape function " << tn << std::endl;
        abort();   
    }
    if (highorder < loworder)
    {
        std::cout << "Error in 'field' object: in 'setorder' the lowest order cannot be larger than the highest order" << std::endl;
        abort();   
    }

    rawfieldptr->setorder(criterion, loworder, highorder, -1); 
}

void field::setorder(double targeterror, int loworder, int highorder, double absthres)
{
    errorifpointerisnull();

    std::string tn = rawfieldptr->gettypename();
    
    if (loworder < 0)
    {
        std::cout << "Error in 'field' object: in 'setorder' cannot use negative minimum interpolation order " << loworder << std::endl;
        abort();   
    }
    hierarchicalformfunction hff;
    if (loworder < hff.getminorder(tn))
    {
        std::cout << "Error in 'field' object: in 'setorder' cannot use interpolation order " << loworder << " for shape function " << tn << std::endl;
        abort();   
    }
    if (highorder < loworder)
    {
        std::cout << "Error in 'field' object: in 'setorder' the lowest order cannot be larger than the highest order" << std::endl;
        abort();   
    }
    
    expression crit = sl::fieldorder(field(getpointer()), 1.0-targeterror, absthres) + 1.0;
    
    crit = sl::max(crit - loworder, 0) + 0.5;
    
    rawfieldptr->setorder(crit, loworder, highorder, highorder-loworder+1);
}

void field::setport(int physreg, port primal, port dual)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
 
    if (primal.getpointer() == dual.getpointer())
    {
        std::cout << "Error in 'field' object: cannot use the same port for the primal and the dual" << std::endl;
        abort();
    }
    
    rawfieldptr->setport(physreg, primal.getpointer(), dual.getpointer());
}

void field::setvalue(int physreg, expression input, int extraintegrationdegree)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setvalue(physreg, -1, NULL, input, extraintegrationdegree);
}

void field::setvalue(int physreg, expression meshdeform, expression input, int extraintegrationdegree)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setvalue(physreg, -1, &meshdeform, input, extraintegrationdegree);
}

void field::setvalue(int physreg, int numfftharms, expression input, int extraintegrationdegree)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setvalue(physreg, numfftharms, NULL, input, extraintegrationdegree);
}

void field::setvalue(int physreg, int numfftharms, expression meshdeform, expression input, int extraintegrationdegree)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setvalue(physreg, numfftharms, &meshdeform, input, extraintegrationdegree);
}

void field::setvalue(int physreg)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setvalue(physreg);
}

void field::setnodalvalues(intdensematrix nodenumbers, densematrix values)
{
    errorifpointerisnull();
    
    if (nodenumbers.count() != values.count())
    {
        std::cout << "Error in 'field' object: argument size mismatch in 'setnodalvalues'" << std::endl;
        abort();
    }
    
    if (nodenumbers.count() > 0)
    {
        int numnodes = universe::mymesh->getelements()->count(0);
        std::vector<int> minmax = nodenumbers.minmax();
        if (minmax[0] < 0)
        {
            std::cout << "Error in 'field' object: trying to set value of node number " << minmax[0] << " (must be positive)" << std::endl;
            abort();
        }
        if (minmax[1] >= numnodes)
        {
            std::cout << "Error in 'field' object: trying to set value of node number " << minmax[1] << " (highest node number in mesh is " << numnodes-1 << ")" << std::endl;
            abort();
        }
    }
    
    rawfieldptr->setnodalvalues(nodenumbers, values);
}

densematrix field::getnodalvalues(intdensematrix nodenumbers)
{
    errorifpointerisnull();
    
    if (nodenumbers.count() > 0)
    {
        int numnodes = universe::mymesh->getelements()->count(0);
        std::vector<int> minmax = nodenumbers.minmax();
        if (minmax[0] < 0)
        {
            std::cout << "Error in 'field' object: trying to get value of node number " << minmax[0] << " (must be positive)" << std::endl;
            abort();
        }
        if (minmax[1] >= numnodes)
        {
            std::cout << "Error in 'field' object: trying to get value of node number " << minmax[1] << " (highest node number in mesh is " << numnodes-1 << ")" << std::endl;
            abort();
        }
    }
    
    return rawfieldptr->getnodalvalues(nodenumbers);
}

void field::setconstraint(int physreg, expression input, int extraintegrationdegree)
{
    errorifpointerisnull(); 
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setdisjregconstraint(physreg, -1, NULL, input, extraintegrationdegree);
}

void field::setconstraint(int physreg, expression meshdeform, expression input, int extraintegrationdegree)
{
    errorifpointerisnull(); 
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setdisjregconstraint(physreg, -1, &meshdeform, input, extraintegrationdegree);
}

void field::setconstraint(int physreg, int numfftharms, expression input, int extraintegrationdegree)
{
    errorifpointerisnull(); 
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setdisjregconstraint(physreg, numfftharms, NULL, input, extraintegrationdegree);
}

void field::setconstraint(int physreg, int numfftharms, expression meshdeform, expression input, int extraintegrationdegree)
{
    errorifpointerisnull(); 
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setdisjregconstraint(physreg, numfftharms, &meshdeform, input, extraintegrationdegree);
}

void field::setconstraint(int physreg)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setdisjregconstraint(physreg);
}

void field::setconditionalconstraint(int physreg, expression condexpr, expression valexpr)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    rawfieldptr->setconditionalconstraint(physreg, condexpr, valexpr);
}

void field::setgauge(int physreg) 
{ 
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    
    if (rawfieldptr->gettypename() != "hcurl")
    {
        std::cout << "Error in 'field' object: cannot gauge shape function " << rawfieldptr->gettypename() << " (only hcurl)" << std::endl;
        abort();   
    }

    rawfieldptr->setgauge(physreg);
}

void field::setdata(int physreg, vectorfieldselect myvec, std::string op) 
{ 
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});
    
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

void field::setcohomologysources(std::vector<int> cutphysregs, std::vector<double> cutvalues)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined(cutphysregs);

    if (cutphysregs.size() != cutvalues.size())
    {
        std::cout << "Error in 'field' object: provided " << cutvalues.size() << " values for " << cutphysregs.size() << " cohomology regions" << std::endl;
        abort();
    }
    
    if (rawfieldptr->gettypename() != "hcurl")
    {
        std::cout << "Error in 'field' object: cannot set a cohomology source to '" << rawfieldptr->gettypename() << "' type fields (use 'hcurl')" << std::endl;
        abort();
    }
    for (int i = 0; i < cutphysregs.size(); i++)
    {
        int prdim = universe::mymesh->getphysicalregions()->get(cutphysregs[i])->getelementdimension();
        if (prdim != -1 && prdim != 1) // -1 for empty is ok
        {
            std::cout << "Error in 'field' object: expected 1D cohomology regions" << std::endl;
            abort();
        }
    }

    if (cutphysregs.size() > 0)
        rawfieldptr->setcohomologysources(cutphysregs, cutvalues);
}

void field::automaticupdate(bool updateit)
{
    errorifpointerisnull();
    rawfieldptr->allowvaluesynchronizing(updateit);
}

void field::noautomaticupdate(void)
{
    errorifpointerisnull();
    rawfieldptr->allowvaluesynchronizing(false);
}

void field::setupdateaccuracy(int extraintegrationorder)
{
    errorifpointerisnull();
    rawfieldptr->setupdateaccuracy(extraintegrationorder);
}

field field::comp(int component) 
{ 
    errorifpointerisnull();
    
    if (component < 0 || component > 2)
    {
        std::cout << "Error in 'field' object: cannot use component number " << component << " (only 0, 1 and 2 are allowed)" << std::endl;
        abort();
    }
    return field(rawfieldptr->comp(component)); 
}

field field::harmonic(const std::vector<int> harmonicnumbers)
{
    errorifpointerisnull();
    
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

void field::writeraw(int physreg, std::string filename, bool isbinary, std::vector<double> extradata)
{
    errorifpointerisnull();
    universe::mymesh->getphysicalregions()->errorundefined({physreg});

    if (isbinary == true && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz" || filename.size() >= 8 && filename.substr(filename.size()-7,7) == ".slz.gz"))
    {
        rawfieldptr->writeraw(physreg, filename, isbinary, extradata);
        return;
    }
    
    if (isbinary == false && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz"))
    {
        rawfieldptr->writeraw(physreg, filename, isbinary, extradata);
        return;
    }
    
    std::cout << "Error in 'field' object: expected .slz file extension (or .slz.gz for binary format) to write raw field data" << std::endl;
    abort();
}

std::vector<double> field::loadraw(std::string filename, bool isbinary)
{
    errorifpointerisnull();

    if (isbinary == true && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz" || filename.size() >= 8 && filename.substr(filename.size()-7,7) == ".slz.gz"))
        return rawfieldptr->loadraw(filename, isbinary);
    
    if (isbinary == false && (filename.size() >= 5 && filename.substr(filename.size()-4,4) == ".slz"))
        return rawfieldptr->loadraw(filename, isbinary);
    
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

expression field::operator+(parameter param) { return (expression)*this + param; }
expression field::operator-(parameter param) { return (expression)*this - param; }
expression field::operator*(parameter param) { return (expression)*this * param; }
expression field::operator/(parameter param) { return (expression)*this / param; }       
 

expression operator+(double val, field inputfield) { return (expression)val + inputfield; }
expression operator-(double val, field inputfield) { return (expression)val - inputfield; }
expression operator*(double val, field inputfield) { return (expression)val * inputfield; }
expression operator/(double val, field inputfield) { return (expression)val / inputfield; }

expression operator+(parameter param, field inputfield) { return (expression)param + inputfield; }
expression operator-(parameter param, field inputfield) { return (expression)param - inputfield; }
expression operator*(parameter param, field inputfield) { return (expression)param * inputfield; }
expression operator/(parameter param, field inputfield) { return (expression)param / inputfield; }

