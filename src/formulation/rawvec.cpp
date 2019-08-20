#include "rawvec.h"


rawvec::rawvec(shared_ptr<dofmanager> dofmngr)
{
    mydofmanager = dofmngr;
    
    VecCreate(PETSC_COMM_SELF, &myvec);
    VecSetSizes(myvec, PETSC_DECIDE, mydofmanager->countdofs());
    VecSetFromOptions(myvec);    
}

rawvec::~rawvec(void) { VecDestroy(&myvec); }

int rawvec::size(void) 
{ 
    if (mydofmanager == NULL)
        return 0;
    else
        return mydofmanager->countdofs(); 
}

void rawvec::removeconstraints(void)
{
    int numdofs = mydofmanager->countdofs();
    int numconstraineddofs = mydofmanager->countconstraineddofs();
    if (numconstraineddofs == 0)
        return;
    
    // Remove the constrained dofs from the dof manager and get the dof renumbering:
    int* dofrenumbering = new int[numdofs];
    mydofmanager = mydofmanager->removeconstraints(dofrenumbering);

    // Rebuild the petsc vector:
    intdensematrix oldaddresses(mydofmanager->countdofs(),1);
    intdensematrix newaddresses(mydofmanager->countdofs(),1);
    int* oldads = oldaddresses.getvalues();
    int* newads = newaddresses.getvalues();
    
    int index = 0;
    for (int i = 0; i < numdofs; i++)
    {
        if (dofrenumbering[i] >= 0)
        {
            oldads[index] = i;
            newads[index] = dofrenumbering[i];
            index++;
        }
    }
    
    densematrix vals = getvalues(oldaddresses);
    
    // Destroy the vector since it will be replaced:
    VecDestroy(&myvec);
    Vec newvec; myvec = newvec;
    
    VecCreate(PETSC_COMM_SELF, &myvec);
    VecSetSizes(myvec, PETSC_DECIDE, mydofmanager->countdofs());
    VecSetFromOptions(myvec);   
    
    setvalues(newaddresses, vals, "set");
    
    delete[] dofrenumbering;
}

void rawvec::updateconstraints(shared_ptr<rawfield> constrainedfield, std::vector<int> disjregs)
{    
    mydofmanager->selectfield(constrainedfield);
    std::vector<shared_ptr<integration>> fieldconstraints = constrainedfield->getconstraints();

    // Loop on all disjoint regions:
    for (int d = 0; d < disjregs.size(); d++)
    {
        int disjreg = disjregs[d];

        // Make sure the disjoint region is in the dof structure and is constrained:
        if (mydofmanager->isdefined(disjreg, 0) && fieldconstraints[disjreg] != NULL)
        {
            // Compute the constraint projection:
            formulation projectconstraint;
            projectconstraint += *fieldconstraints[disjreg];
            projectconstraint.isconstraintcomputation = true;                

            // Get an all zero vector:
            vec constraintvalvec(projectconstraint);
            // Zero valued constraints need not be computed:
            if (fieldconstraints[disjreg]->isprojectionofzero == false)
            {
                projectconstraint.generate();
                constraintvalvec = mathop::solve(projectconstraint.A(), projectconstraint.b());
            }

            // Loop on all disjoint regions who share the same constraint-computation-formulation:
            shared_ptr<integration> currentconstraint = fieldconstraints[disjreg];
            for (int i = d; i < disjregs.size(); i++)
            {
                if (mydofmanager->isdefined(disjregs[i], 0) && fieldconstraints[disjregs[i]] == currentconstraint)
                {
                    // Transfer the data from 'constraintvalvec' to 'myvec' on disjoint region disjregs[i]:
                    setdata(constraintvalvec.getpointer(), disjregs[i], constrainedfield);
                    // Clear the fieldconstraints pointer since it was treated:
                    fieldconstraints[disjregs[i]] = NULL;
                }
            }
        }
    }
}

void rawvec::setvalues(intdensematrix addresses, densematrix valsmat, std::string op)
{            
    double* myval = valsmat.getvalues();
    int* myad = addresses.getvalues();

    long int numentries = addresses.count();

    // Remove the negative addresses:
    long int numpositiveentries = addresses.countpositive();
    
    if (numpositiveentries == 0)
    	return;

    intdensematrix filteredad(numpositiveentries,1);
    densematrix filteredval(numpositiveentries,1);
    int* filteredads = filteredad.getvalues();
    double* filteredvals = filteredval.getvalues();

    long int index = 0;
    for (long int i = 0; i < numentries; i++)
    {
        if (myad[i] >= 0)
        {
            filteredads[index] = myad[i];
            filteredvals[index] = myval[i];
            index++;        
        }
    }

    if (op == "add")
        VecSetValues(myvec, numpositiveentries, filteredads, filteredvals, ADD_VALUES);
    if (op == "set")
        VecSetValues(myvec, numpositiveentries, filteredads, filteredvals, INSERT_VALUES);
}

void rawvec::setvalue(int address, double value, std::string op)
{            
    if (address < 0)
        return;
    
    if (op == "add")
        VecSetValue(myvec, address, value, ADD_VALUES);
    if (op == "set")
        VecSetValue(myvec, address, value, INSERT_VALUES);
}

densematrix rawvec::getvalues(intdensematrix addresses)
{
    long int numentries = addresses.count();
    densematrix valmat(numentries,1);
    VecGetValues(myvec, numentries, addresses.getvalues(), valmat.getvalues());
    
    return valmat;
}

double rawvec::getvalue(int address)
{
	int ads[1] = {address};
	double outval[1];
	VecGetValues(myvec, 1, ads, outval);
	
	return outval[0];
}

void rawvec::setvalues(shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex, densematrix vals, std::string op)
{
    mydofmanager->selectfield(selectedfield);

    // Do nothing if the entries are not in the vector:
    if (mydofmanager->isdefined(disjointregionnumber, formfunctionindex) == false)
        return;

    int rangebegin = mydofmanager->getrangebegin(disjointregionnumber, formfunctionindex);
    int rangeend = mydofmanager->getrangeend(disjointregionnumber, formfunctionindex);

    int numentries = rangeend-rangebegin+1;

    intdensematrix addressestoset(numentries, 1, rangebegin, 1);
    setvalues(addressestoset, vals, op);
}

densematrix rawvec::getvalues(shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex)
{
    mydofmanager->selectfield(selectedfield);
    
    // Return an empty matrix if the entries are not in the vector:
    if (mydofmanager->isdefined(disjointregionnumber, formfunctionindex) == false)
    {
    	densematrix vals;
    	return vals;
    }
    
    int rangebegin = mydofmanager->getrangebegin(disjointregionnumber, formfunctionindex);
    int rangeend = mydofmanager->getrangeend(disjointregionnumber, formfunctionindex);
    
    int numentries = rangeend-rangebegin+1;
    
    densematrix vals(numentries, 1);
    intdensematrix addressestoget(numentries, 1, rangebegin, 1);
    VecGetValues(myvec, numentries, addressestoget.getvalues(), vals.getvalues());
    
    return vals;
}

void rawvec::write(std::string filename)
{
    if (mydofmanager == NULL)
    {
        std::cout << "Error in 'rawvec' object: cannot write vec (structure is not defined)" << std::endl;
        abort();
    }
    
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);
        
        if (fileext == ".txt")
        {
            std::vector<double> towrite;
            intdensematrix adresses(1,size(),0,1);
            densematrix curvals = getvalues(adresses);
            curvals.getvalues(towrite);
            mathop::writevector(filename, towrite, ',', true);
            return;
        }
        if (fileext == ".bin")
        {
            PetscViewer v;
            PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_WRITE, &v);
            VecView(myvec, v);
            PetscViewerDestroy(&v);
            return;
        }
    }
    
    std::cout << "Error in 'rawvec' object: cannot write vec data to file '" << filename << "'." << std::endl;
    std::cout << "Supported formats are .txt (ASCII) and .bin (binary)." << std::endl;
    abort();
}

void rawvec::load(std::string filename)
{
    if (mydofmanager == NULL)
    {
        std::cout << "Error in 'rawvec' object: cannot load vec (structure is not defined)" << std::endl;
        abort();
    }
    
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);
        
        if (fileext == ".txt")
        {
            std::vector<double> loadedvals = mathop::loadvector(filename, ',', true);
            if (loadedvals.size() != size())
            {
                std::cout << "Error in 'rawvec' object: loaded data size (" << loadedvals.size() << ") does not match the vector size (" << size() << ") " << std::endl;
                abort();
            }
            densematrix valsmat(1,size(), loadedvals);
            intdensematrix adresses(1,size(),0,1);
            setvalues(adresses, valsmat);
            return;
        }
        if (fileext == ".bin")
        {
            PetscViewer v;
            PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_READ, &v);
            VecLoad(myvec, v);
            PetscViewerDestroy(&v);
            return;
        }
    }
    
    std::cout << "Error in 'rawvec' object: cannot load vec data from file '" << filename << "'." << std::endl;
    std::cout << "Supported formats are .txt (ASCII) and .bin (binary)." << std::endl;
    abort();
}
     
void rawvec::print(void)
{          
    std::cout << std::endl;
    if (mydofmanager == NULL)
        std::cout << "Undefined vector" << std::endl;
    else
    {
        if (size() > 0)
            VecView(myvec, PETSC_VIEWER_STDOUT_SELF);
        else
            std::cout << "Empty vector" << std::endl;
    }
    std::cout << std::endl;
}

void rawvec::setdata(shared_ptr<rawvec> inputvec, int disjreg, shared_ptr<rawfield> inputfield)
{
    mydofmanager->selectfield(inputfield);
    inputvec->mydofmanager->selectfield(inputfield);
        
    int ff = 0;
    while (mydofmanager->isdefined(disjreg, ff) && inputvec->mydofmanager->isdefined(disjreg, ff))
    {
        int myrangebegin = mydofmanager->getrangebegin(disjreg,ff);
        int myrangeend = mydofmanager->getrangeend(disjreg,ff);
        
        int inputrangebegin = inputvec->mydofmanager->getrangebegin(disjreg,ff);
        
        int numdofs = myrangeend-myrangebegin+1;
        
        intdensematrix myaddresses(numdofs, 1, myrangebegin, 1);
        intdensematrix inputaddresses(numdofs, 1, inputrangebegin, 1);
        
        densematrix inputval = inputvec->getvalues(inputaddresses);
        setvalues(myaddresses, inputval, "set");
        
        ff++;
    }
    
    // In case there are too few form functions in 'inputvec' set a 0 value.
    while (mydofmanager->isdefined(disjreg, ff))
    {
        int myrangebegin = mydofmanager->getrangebegin(disjreg,ff);
        int myrangeend = mydofmanager->getrangeend(disjreg,ff);
                
        int numdofs = myrangeend-myrangebegin+1;
        
        intdensematrix myaddresses(numdofs, 1, myrangebegin, 1);
        
        densematrix zerovals(numdofs,1, 0);
        setvalues(myaddresses, zerovals, "set");
        
        ff++;
    }
}




