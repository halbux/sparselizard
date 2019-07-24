#include "vec.h"


vec::vec(formulation formul) { rawvecptr = shared_ptr<rawvec>(new rawvec(formul.getdofmanager())); }

void vec::errorifpointerisnull(void)
{
    if (rawvecptr == NULL)
    {
        std::cout << "Error in 'vec' object: cannot perform the operation (vector is undefined)" << std::endl;
        abort();
    }
}

int vec::size(void) { errorifpointerisnull(); return rawvecptr->size(); }

void vec::removeconstraints(void) { errorifpointerisnull(); rawvecptr->removeconstraints(); };
        
void vec::updateconstraints(void)
{
    errorifpointerisnull();
    
    std::vector<int> disjregs((universe::mymesh->getdisjointregions())->count());
	// Set 'disjregs' to [0 1 2 ...]:
   	std::iota(disjregs.begin(), disjregs.end(), 0);
    
    std::vector<shared_ptr<rawfield>> fieldsindofmanager = rawvecptr->getdofmanager()->getfields();
    for (int i = 0; i < fieldsindofmanager.size(); i++)
        rawvecptr->updateconstraints(fieldsindofmanager[i], disjregs);
        
    // Update the conditional constraints:
    shared_ptr<dofmanager> mydofmanager = rawvecptr->getdofmanager();
    std::pair<intdensematrix, densematrix> condconstrdata = mydofmanager->getconditionalconstraintdata();
    rawvecptr->setvalues(condconstrdata.first, condconstrdata.second);
    
    // Set the gauged indexes to zero:
    intdensematrix gaugedindexes = mydofmanager->getgaugedindexes();
    int numgaugedindexes = gaugedindexes.count();
    if (numgaugedindexes > 0)
        rawvecptr->setvalues(gaugedindexes, densematrix(gaugedindexes.countrows(),gaugedindexes.countcolumns(), 0.0));
}

void vec::setvalues(intdensematrix addresses, densematrix valsmat, std::string op) 
{ 
	errorifpointerisnull(); rawvecptr->setvalues(addresses, valsmat, op); 
}
densematrix vec::getvalues(intdensematrix addresses) { errorifpointerisnull(); return rawvecptr->getvalues(addresses); }

void vec::setvalue(int address, double value, std::string op)
{
	errorifpointerisnull(); rawvecptr->setvalue(address, value, op); 
}

double vec::getvalue(int address)
{
	errorifpointerisnull(); return rawvecptr->getvalue(address); 
}

vectorfieldselect vec::operator|(field selectedfield) { errorifpointerisnull(); return vectorfieldselect(rawvecptr, selectedfield.getpointer()); }

void vec::setdata(int physreg, field myfield, std::string op)
{
    errorifpointerisnull();
    vectorfieldselect(rawvecptr, myfield.getpointer()).setdata(physreg, myfield, op);
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
    return vec(shared_ptr<rawvec>(new rawvec(  rawvecptr->getdofmanager(), output  )));
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


vec vec::operator+(void) { return copy(); }
vec vec::operator-(void) { return *this*-1; }

vec vec::operator*(double input)
{
    Vec x = getpetsc();
    Vec output;
    VecDuplicate(x, &output);
    VecAXPY(output, input, x);
    return vec(shared_ptr<rawvec>(new rawvec(  rawvecptr->getdofmanager(), output  )));
}

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
