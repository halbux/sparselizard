#include "rawvec.h"


void rawvec::synchronize(void)
{
    if (issynchronizing || mymeshtracker == universe::mymesh->getmeshtracker() || mydofmanager == NULL)
        return;
    issynchronizing = true; 

    
    // Extract all the values from the vector:
    intdensematrix adresses(mycurrentstructure[0].countdofs(),1, 0,1);
    densematrix alloldvals = getvalues(adresses);
    double* alloldvalsptr = alloldvals.getvalues();
    
    // Create new petsc object (needed since number of dofs can change):
    VecDestroy(&myvec);
    VecCreate(PETSC_COMM_SELF, &myvec);
    VecSetSizes(myvec, PETSC_DECIDE, mydofmanager->countdofs());
    VecSetFromOptions(myvec);    
    

    // We need to know in which disjoint region every old element number was:
    std::vector<std::vector<int>> inolddisjregs;
    mymeshtracker->getindisjointregions(inolddisjregs);
    // We also need to know how the mesh structure was before:
    disjointregions* olddisjregs = mymeshtracker->getdisjointregions();

    
    // Get the map to renumber the current elements to the ones in the mesh as in the current rawvec status:
    std::vector<std::vector<int>> numberback;
    universe::mymesh->getmeshtracker()->getrenumbering(mymeshtracker, numberback);
    
    // The fields are unchanged:
    std::vector<std::shared_ptr<rawfield>> curfields = mydofmanager->getfields();


    for (int i = 0; i < curfields.size(); i++)
    {
        mycurrentstructure[0].selectfield(curfields[i]);
        mydofmanager->selectfield(curfields[i]);
    
        for (int d = 0; d < universe::mymesh->getdisjointregions()->count(); d++)
        {
            int numff = mydofmanager->countformfunctions(d);
            if (numff == 0)
                continue;
            
            int elemtypenum = universe::mymesh->getdisjointregions()->getelementtypenumber(d);
            int rb = universe::mymesh->getdisjointregions()->getrangebegin(d);
        
            for (int ff = 0; ff < numff; ff++)
            {
                densematrix vals(universe::mymesh->getdisjointregions()->countelements(d), 1, 0.0);
                double* myvals = vals.getvalues();
            
                for (int e = 0; e < vals.countrows(); e++)
                {
                    // Find the corresponding index and value in 'alloldvalsptr':
                    int oldelemnum = numberback[elemtypenum][rb+e];
                    int oldelemdisjreg = inolddisjregs[elemtypenum][oldelemnum];
                    int oldelemindex = oldelemnum - olddisjregs->getrangebegin(oldelemdisjreg);
                    
                    // Value will be 0 if not in the previous vector:
                    if (mycurrentstructure[0].isdefined(oldelemdisjreg, ff))
                    {
                        int olddofrb = mycurrentstructure[0].getrangebegin(oldelemdisjreg, ff);
                        myvals[e] = alloldvalsptr[olddofrb+oldelemindex];
                    }
                }
                
                setvalues(curfields[i], d, ff, vals, "set");
            }
        }
    }    
    
    
    // Update the dof manager to the current one:
    mycurrentstructure = {*mydofmanager};
    mycurrentstructure[0].donotsynchronize();
    
    // Update the mesh tracker to the current one:
    mymeshtracker = universe::mymesh->getmeshtracker();
    issynchronizing = false;
}

rawvec::rawvec(std::shared_ptr<dofmanager> dofmngr)
{
    mydofmanager = dofmngr;

    VecCreate(PETSC_COMM_SELF, &myvec);
    VecSetSizes(myvec, PETSC_DECIDE, mydofmanager->countdofs());
    VecSetFromOptions(myvec);   
    
    mymeshtracker = universe::mymesh->getmeshtracker();
    
    mycurrentstructure = {*dofmngr};
    mycurrentstructure[0].donotsynchronize();
}

rawvec::rawvec(std::shared_ptr<dofmanager> dofmngr, Vec input)
{
    mydofmanager = dofmngr;

    myvec = input;
    
    mymeshtracker = universe::mymesh->getmeshtracker();
    
    mycurrentstructure = {*dofmngr};
    mycurrentstructure[0].donotsynchronize();
}

rawvec::~rawvec(void) { VecDestroy(&myvec); }

int rawvec::size(void) 
{ 
    synchronize();

    if (mydofmanager == NULL)
        return 0;
    else
        return mydofmanager->countdofs(); 
}

void rawvec::removeconstraints(void)
{
    synchronize();
    
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
    
    // The current structure tracker must be updated to the new dof manager:
    mycurrentstructure = {*mydofmanager};
    mycurrentstructure[0].donotsynchronize();
}

void rawvec::updateconstraints(std::shared_ptr<rawfield> constrainedfield, std::vector<int> disjregs)
{    
    synchronize();
    
    mydofmanager->selectfield(constrainedfield);
    std::vector<std::shared_ptr<integration>> fieldconstraints = constrainedfield->getconstraints();

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
            std::shared_ptr<integration> currentconstraint = fieldconstraints[disjreg];
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
    synchronize();
     
    double* myval = valsmat.getvalues();
    int* myad = addresses.getvalues();

    int numentries = addresses.count();

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
    synchronize();
    
    if (address < 0)
        return;
    
    if (op == "add")
        VecSetValue(myvec, address, value, ADD_VALUES);
    if (op == "set")
        VecSetValue(myvec, address, value, INSERT_VALUES);
}

densematrix rawvec::getvalues(intdensematrix addresses)
{
    synchronize();
    
    int numentries = addresses.count();
    densematrix valmat(numentries,1);
    VecGetValues(myvec, numentries, addresses.getvalues(), valmat.getvalues());
    
    return valmat;
}

double rawvec::getvalue(int address)
{
    synchronize();
    
    int ads[1] = {address};
    double outval[1];
    VecGetValues(myvec, 1, ads, outval);
    
    return outval[0];
}

void rawvec::setvalues(std::shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex, densematrix vals, std::string op)
{
    synchronize();
    
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

densematrix rawvec::getvalues(std::shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex)
{
    synchronize();
    
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
    synchronize();
    
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
    synchronize();
    
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
    synchronize();
     
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

std::shared_ptr<dofmanager> rawvec::getdofmanager(void)
{
    synchronize();
    
    return mydofmanager;
}

Vec rawvec::getpetsc(void)
{
    synchronize();
    
    return myvec;
}

void rawvec::setdata(std::shared_ptr<rawvec> inputvec, int disjreg, std::shared_ptr<rawfield> inputfield)
{
    synchronize();
    
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

