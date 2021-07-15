#include "vec.h"


vec::vec(formulation formul) { rawvecptr = std::shared_ptr<rawvec>(new rawvec(formul.getdofmanager())); }

void vec::errorifpointerisnull(void)
{
    if (rawvecptr == NULL)
    {
        std::cout << "Error in 'vec' object: cannot perform the operation (vector is undefined)" << std::endl;
        abort();
    }
}

vec::vec(int vecsize, intdensematrix addresses, densematrix vals)
{
    rawvecptr = std::shared_ptr<rawvec>(new rawvec(std::shared_ptr<dofmanager>(new dofmanager(vecsize))));
    rawvecptr->setvalues(addresses, vals, "set");
}

int vec::size(void) { errorifpointerisnull(); return rawvecptr->size(); }

void vec::permute(intdensematrix rowpermute, bool invertit)
{
    if (rowpermute.count() != size())
    {
        std::cout << "Error in 'vec' object: unexpected argument size for permutation" << std::endl;
        abort();
    }

    IS rowpermutis;
    ISCreateGeneral(PETSC_COMM_SELF, rowpermute.count(), rowpermute.getvalues(), PETSC_USE_POINTER, &rowpermutis);
    ISSetPermutation(rowpermutis);
    
    if (invertit == false)
        VecPermute(getpetsc(), rowpermutis, PETSC_FALSE);
    else
        VecPermute(getpetsc(), rowpermutis, PETSC_TRUE);
}

void vec::updateconstraints(void)
{
    errorifpointerisnull();
    
    std::shared_ptr<dofmanager> mydofmanager = rawvecptr->getdofmanager();
    // Get all disjoint regions:
    std::vector<int> disjregs((universe::mymesh->getdisjointregions())->count());
    std::iota(disjregs.begin(), disjregs.end(), 0);
    
    // Update the disjoint region constraints:
    std::vector<std::shared_ptr<rawfield>> fieldsindofmanager = rawvecptr->getdofmanager()->getfields();
    for (int i = 0; i < fieldsindofmanager.size(); i++)
        rawvecptr->updatedisjregconstraints(fieldsindofmanager[i], disjregs);
        
    // Update the conditional constraints:
    std::pair<intdensematrix, densematrix> condconstrdata = mydofmanager->getconditionalconstraintdata();
    rawvecptr->setvalues(condconstrdata.first, condconstrdata.second);
    
    // Set the gauged indexes to zero:
    intdensematrix gaugedindexes = mydofmanager->getgaugedindexes();
    rawvecptr->setvalues(gaugedindexes, densematrix(gaugedindexes.countrows(), gaugedindexes.countcolumns(), 0.0));
}

void vec::setvalues(intdensematrix addresses, densematrix valsmat, std::string op) 
{ 
    errorifpointerisnull(); rawvecptr->setvalues(addresses, valsmat, op); 
}

void vec::setallvalues(densematrix valsmat, std::string op)
{ 
    errorifpointerisnull();
    intdensematrix ads(size(),1,0,1);
    rawvecptr->setvalues(ads, valsmat, op); 
}

densematrix vec::getvalues(intdensematrix addresses) { errorifpointerisnull(); return rawvecptr->getvalues(addresses); }

densematrix vec::getallvalues(void)
{
    errorifpointerisnull();
    intdensematrix ads(size(),1,0,1);
    return rawvecptr->getvalues(ads);
}

void vec::setvalue(int address, double value, std::string op)
{
    errorifpointerisnull(); rawvecptr->setvalue(address, value, op); 
}

double vec::getvalue(int address)
{
    errorifpointerisnull(); return rawvecptr->getvalue(address); 
}

void vec::setvalue(port prt, double value, std::string op)
{
    errorifpointerisnull();
    
    if (prt.getpointer()->isharmonicone())
    {
        int address = rawvecptr->getdofmanager()->getaddress(prt.getpointer()->harmonic(1).get());
        setvalue(address, value, op);
    }
    else
    {
        std::cout << "Error in 'vec' object: cannot set the value of a multiharmonic port (only constant harmonic 1)" << std::endl;
        abort();
    }
}

double vec::getvalue(port prt)
{
    errorifpointerisnull();
    
    if (prt.getpointer()->isharmonicone())
    {
        int address = rawvecptr->getdofmanager()->getaddress(prt.getpointer()->harmonic(1).get());
        return getvalue(address);
    }
    else
    {
        std::cout << "Error in 'vec' object: cannot get the value of a multiharmonic port (only constant harmonic 1)" << std::endl;
        abort();
    }
}

vectorfieldselect vec::operator|(field selectedfield) { errorifpointerisnull(); return vectorfieldselect(rawvecptr, selectedfield.getpointer()); }

void vec::setdata(int physreg, field myfield, std::string op)
{
    errorifpointerisnull();
    vectorfieldselect(rawvecptr, myfield.getpointer()).setdata(physreg, myfield, op);
}

void vec::setdata(void)
{
    errorifpointerisnull();
    
    std::vector<std::shared_ptr<rawfield>> rfs = rawvecptr->getdofmanager()->getfields();
    for (int i = 0; i < rfs.size(); i++)
        rfs[i]->transferdata(-1, vectorfieldselect(rawvecptr, rfs[i]), "set");
        
    rawvecptr->setvaluesfromports();
}

void vec::automaticupdate(bool updateit)
{
    errorifpointerisnull();
    rawvecptr->allowvaluesynchronizing(updateit);
}

void vec::noautomaticupdate(void)
{
    errorifpointerisnull();
    rawvecptr->allowvaluesynchronizing(false);
}

Vec vec::getpetsc(void) { errorifpointerisnull(); return rawvecptr->getpetsc(); }

void vec::write(std::string filename)
{
    errorifpointerisnull();
    rawvecptr->write(filename);
}

void vec::load(std::string filename)
{
    errorifpointerisnull();
    rawvecptr->load(filename);
}

void vec::print(void) { errorifpointerisnull(); rawvecptr->print(); }

vec vec::copy(void)
{
    Vec x = getpetsc();
    Vec output;
    VecDuplicate(x, &output);
    VecCopy(x, output);
    return vec(std::shared_ptr<rawvec>(new rawvec(  rawvecptr->getdofmanager(), output  )));
}

vec vec::extract(intdensematrix addresses)
{
    densematrix extractedvals = getvalues(addresses);
    std::shared_ptr<rawvec> newrawvecptr(new rawvec(std::shared_ptr<dofmanager>(new dofmanager(addresses.count()))));
    newrawvecptr->setvalues(intdensematrix(addresses.count(), 1, 0, 1), extractedvals, "set");
    return vec(newrawvecptr);
}

double vec::norm(std::string type)
{
    double normval;
    
    if (type == "1") { VecNorm(getpetsc(), NORM_1, &normval); return normval; }
    if (type == "2") { VecNorm(getpetsc(), NORM_2, &normval); return normval; }
    if (type == "infinity") { VecNorm(getpetsc(), NORM_INFINITY, &normval); return normval; }
    
    std::cout << "Error in 'vec' object: norm type unknown (use '1', '2' or 'infinity)" << std::endl;
    abort();
}

double vec::sum(void)
{
    double sumval;
    VecSum(getpetsc(), &sumval);
    return sumval;
}


vec vec::operator+(void) { return copy(); }
vec vec::operator-(void) { return *this*-1; }

vec vec::operator*(double input)
{
    Vec x = getpetsc();
    Vec output;
    VecDuplicate(x, &output);
    VecAXPY(output, input, x);
    return vec(std::shared_ptr<rawvec>(new rawvec(  rawvecptr->getdofmanager(), output  )));
}

vec vec::operator/(double input) { return *this*(1.0/input); }

vec vec::operator+(vec input)
{
    vec copied = copy();
    Vec output = copied.getpetsc();
    Vec x = input.getpetsc();
    VecAXPY(output, 1, x);
    return copied;
}

vec vec::operator-(vec input)
{
    vec copied = copy();
    Vec output = copied.getpetsc();
    Vec x = input.getpetsc();
    VecAXPY(output, -1, x);
    return copied;
}


vec operator*(double inputdouble, vec inputvec) { return inputvec*inputdouble; }

