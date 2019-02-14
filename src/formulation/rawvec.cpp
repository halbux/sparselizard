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

    int numentries = addresses.countrows()*addresses.countcolumns();

    // Remove the negative addresses:
    int numpositiveentries = addresses.countpositive();
    
    if (numpositiveentries == 0)
    	return;

    intdensematrix filteredad(numpositiveentries,1);
    densematrix filteredval(numpositiveentries,1);
    int* filteredads = filteredad.getvalues();
    double* filteredvals = filteredval.getvalues();

    int index = 0;
    for (int i = 0; i < numentries; i++)
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
    int numentries = addresses.countrows()*addresses.countcolumns();
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




