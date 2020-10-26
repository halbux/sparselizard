#include "rawvec.h"


void rawvec::synchronize(void)
{
    if (mydofmanager == NULL || mydofmanager->ismanaged() == false || issynchronizing || myptracker == universe::mymesh->getptracker())
        return;
    issynchronizing = true; 

    std::vector<std::shared_ptr<rawfield>> dmfields = mycurrentstructure[0].getfields();
    std::vector<std::shared_ptr<rawfield>> datafields(dmfields.size());
    
    if (isvaluesynchronizingallowed)
    {
        // For a correct 'setdata' call below (now 'mydofmanager' will not be synced either):
        dofmanager dm;
        dm = *mydofmanager; // backup
        *mydofmanager = mycurrentstructure[0];
        
        for (int i = 0; i < dmfields.size(); i++)
        {
            mydofmanager->selectfield(dmfields[i]);
            mycurrentstructure[0].selectfield(dmfields[i]);
            datafields[i] = std::shared_ptr<rawfield>(new rawfield(&(mycurrentstructure[0]), myrawmesh, myptracker));
            mydofmanager->replaceselectedfield(datafields[i]);
            mycurrentstructure[0].replaceselectedfield(datafields[i]);
            datafields[i]->allowsynchronizing(false);
            datafields[i]->setdata(-1, vec(shared_from_this())|field(datafields[i]));
            datafields[i]->allowsynchronizing(true);
        }
          
        *mydofmanager = dm; // restore
        
        for (int i = 0; i < dmfields.size(); i++)
        {
            mydofmanager->selectfield(dmfields[i]);
            datafields[i]->synchronize({}, mydofmanager->getselectedfieldorders());
        }
    }
    
    // Create new petsc object (needed since number of dofs can change):
    VecDestroy(&myvec);
    VecCreate(PETSC_COMM_SELF, &myvec);
    VecSetSizes(myvec, PETSC_DECIDE, mydofmanager->countdofs());
    VecSetFromOptions(myvec);    
    
    // Update the dof manager to the current one:
    mycurrentstructure = {*mydofmanager};
    mycurrentstructure[0].donotsynchronize();
    
    // Update the mesh tracker to the current one:
    myptracker = universe::mymesh->getptracker();
    myrawmesh = universe::mymesh;
    
    if (isvaluesynchronizingallowed)
    {
        // Transfer the data back to the vector:
        for (int i = 0; i < datafields.size(); i++)
            datafields[i]->transferdata(-1, vec(shared_from_this())|field(dmfields[i]), "set");
    }
    
    issynchronizing = false;
}

rawvec::rawvec(std::shared_ptr<dofmanager> dofmngr)
{
    mydofmanager = dofmngr;

    VecCreate(PETSC_COMM_SELF, &myvec);
    VecSetSizes(myvec, PETSC_DECIDE, mydofmanager->countdofs());
    VecSetFromOptions(myvec);   
    
    if (mydofmanager->ismanaged())
    {
        myptracker = universe::mymesh->getptracker();
        myrawmesh = universe::mymesh;
        
        mycurrentstructure = {*dofmngr};
        mycurrentstructure[0].donotsynchronize();
    }
}

rawvec::rawvec(std::shared_ptr<dofmanager> dofmngr, Vec input)
{
    mydofmanager = dofmngr;

    myvec = input;
    
    if (mydofmanager->ismanaged())
    {
        myptracker = universe::mymesh->getptracker();
        myrawmesh = universe::mymesh;
        
        mycurrentstructure = {*dofmngr};
        mycurrentstructure[0].donotsynchronize();
    }
}

rawvec::~rawvec(void)
{
    // Avoid crashes when destroy is called after PetscFinalize (not allowed).
    PetscBool ispetscinitialized;
    PetscInitialized(&ispetscinitialized);

    if (ispetscinitialized == PETSC_TRUE)
        VecDestroy(&myvec);
}

void rawvec::allowvaluesynchronizing(bool allowit)
{
    isvaluesynchronizingallowed = allowit;
}

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

std::shared_ptr<rawvec> rawvec::getpointer(void)
{
    synchronize();
   
    return shared_from_this();
}

std::shared_ptr<rawmesh> rawvec::getrawmesh(void)
{
    synchronize();
    
    return myrawmesh;
}

std::shared_ptr<ptracker> rawvec::getptracker(void)
{
    synchronize();
    
    return myptracker;
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

