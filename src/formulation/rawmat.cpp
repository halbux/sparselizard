#include "rawmat.h"


rawmat::~rawmat(void) 
{ 
    MatDestroy(&mymat);
    if (ludefined) 
        KSPDestroy(&myksp);
}

int rawmat::countrows(void) 
{ 
    if (mydofmanager == NULL)
        return 0;
    else
        return mydofmanager->countdofs();
}

int rawmat::countcolumns(void) 
{ 
    if (mydofmanager == NULL)
        return 0;
    else
        return mydofmanager->countdofs();
}

void rawmat::zeroentries(intdensematrix entriestozero, bool zerorows, bool zerocolumns)
{
    int* entriestozeroptr = entriestozero.getvalues();

    int numentries = entriestozero.count();
    if (numentries > 0)
    {
        // First build a vector for direct access:
        std::vector<bool> istozero(mydofmanager->countdofs(), false);

        for (int i = 0; i < numentries; i++)
            istozero[entriestozeroptr[i]] = true;
            
        if (zerorows)
        {
            for (int i = 0; i < accumulatedrowindices.size(); i++)
            {
                int* accumulatedrowindicesptr = accumulatedrowindices[i].getvalues();
                for (int j = 0; j < accumulatedrowindices[i].count(); j++)
                {
                    if (accumulatedrowindicesptr[j] >= 0 && istozero[accumulatedrowindicesptr[j]])
                        accumulatedrowindicesptr[j] = -1;
                }
            }
        }
        if (zerocolumns)
        {
            for (int i = 0; i < accumulatedcolindices.size(); i++)
            {
                int* accumulatedcolindicesptr = accumulatedcolindices[i].getvalues();
                for (int j = 0; j < accumulatedcolindices[i].count(); j++)
                {
                    if (accumulatedcolindicesptr[j] >= 0 && istozero[accumulatedcolindicesptr[j]])
                        accumulatedcolindicesptr[j] = -1;
                }
            }
        }
    }
}

void rawmat::gauge(void)
{
    intdensematrix gaugedindexes = mydofmanager->getgaugedindexes();
    zeroentries(gaugedindexes, true, true);
}


void rawmat::removeconstraints(void)
{
    if (mydofmanager->countconstraineddofs() == 0)
        return;
    
    if (ludefined)
    {
        std::cout << "Error in 'rawmat' object: cannot remove the constraints when an LU decomposition is associated" << std::endl;
        abort();
    }
    if (petscrows.getvalues() == NULL)
    {
        std::cout << "Error in 'rawmat' object: cannot remove the constraints for this matrix (make sure you got the matrix directly from the formulation and it is not a copy)" << std::endl;
        abort();
    }
    
    // Free the matrix since it will be replaced:
    if (nnz >= 0) 
    {
        MatDestroy(&mymat);
        Mat newmat; mymat = newmat;
    }
        
    // Bring 'petscrows' from CSR to ijk format:
    intdensematrix ijkpetscrows(petscvals.count(),1);
    myalgorithm::csrtoijk(mydofmanager->countdofs(), petscrows.getvalues(), ijkpetscrows.getvalues());
    petscrows = ijkpetscrows;
    
    // Remove the constrained dofs from the dof manager and get the dof renumbering:
    int* dofrenumbering = new int[mydofmanager->countdofs()];
    mydofmanager = mydofmanager->removeconstraints(dofrenumbering);

    
    // Rebuild the petsc matrix:
    
    accumulatedrowindices = {petscrows};
    accumulatedcolindices = {petsccols};
    accumulatedvals = {petscvals};
    
    int* myrows = petscrows.getvalues();
    int* mycols = petsccols.getvalues();
    
    for (int i = 0; i < petscrows.count(); i++)
    {
        myrows[i] = dofrenumbering[myrows[i]];
        mycols[i] = dofrenumbering[mycols[i]];
    }
    
    delete[] dofrenumbering;
    
    process();
    clearfragments();
}

void rawmat::accumulate(intdensematrix rowadresses, intdensematrix coladresses, densematrix vals)
{
    accumulatedrowindices.push_back(rowadresses);
    accumulatedcolindices.push_back(coladresses);
    accumulatedvals.push_back(vals);
}

void rawmat::process(void)
{
    // Concatenate all accumulated fragments!
    // Remove negative indexes. First get the overall length.
    int veclen = 0;
    for (int i = 0; i < accumulatedvals.size(); i++)
    {
        int* accumulatedrowindicesptr = accumulatedrowindices[i].getvalues();
        int* accumulatedcolindicesptr = accumulatedcolindices[i].getvalues();
        
        for (int j = 0; j < accumulatedvals[i].count(); j++)
        {
            if (accumulatedrowindicesptr[j] >= 0 && accumulatedcolindicesptr[j] >= 0)
                veclen++;            
        }
    }
    
    // Stitch all fragments together while removing negative indexes.
    int* stitchedrowindices = new int[veclen];
    int* stitchedcolindices = new int[veclen];
    double* stitchedvals = new double[veclen];

    int ind = 0;
    for (int i = 0; i < accumulatedvals.size(); i++)
    {
        double* valsptr = accumulatedvals[i].getvalues();
        int* accumulatedrowindicesptr = accumulatedrowindices[i].getvalues();
        int* accumulatedcolindicesptr = accumulatedcolindices[i].getvalues();
        
        for (int j = 0; j < accumulatedvals[i].count(); j++)
        {
            if (accumulatedrowindicesptr[j] >= 0 && accumulatedcolindicesptr[j] >= 0)
            {
                stitchedrowindices[ind] = accumulatedrowindicesptr[j];
                stitchedcolindices[ind] = accumulatedcolindicesptr[j];
                stitchedvals[ind] = valsptr[j];
                ind++;            
            }
        }
    }
    
    // Sort the vectors according to the row then according to the column. 
    // 'reorderingvector' will tell us how to reorder.
    int* reorderingvector = myalgorithm::stablesortparallel({stitchedrowindices, stitchedcolindices}, veclen);
        
    // Get the number of nonzeros:
    nnz = 0;
    if (veclen > 0)
        nnz = 1;
    for (int i = 1; i < veclen; i++)
    {
        if (stitchedrowindices[reorderingvector[i]] != stitchedrowindices[reorderingvector[i-1]] || stitchedcolindices[reorderingvector[i]] != stitchedcolindices[reorderingvector[i-1]])
            nnz++;
    }

    // We now add together the values at the same adresses. Sparse format used is CSR.
    petscrows = intdensematrix(countrows()+1,1);
    petsccols = intdensematrix(nnz,1);
    petscvals = densematrix(nnz,1);
    
    int* finalrowindices = petscrows.getvalues();
    int* finalcolindices = petsccols.getvalues();
    double* finalvals = petscvals.getvalues();

    int row = -1; 
    ind = -1;
    for (int i = 0; i < veclen; i++)
    {
        // Same row:
        if (i > 0 && stitchedrowindices[reorderingvector[i]] == stitchedrowindices[reorderingvector[i-1]])
        {
            // Same row and column:
            if (stitchedcolindices[reorderingvector[i]] == stitchedcolindices[reorderingvector[i-1]])
                finalvals[ind] += stitchedvals[reorderingvector[i]];
            else
            {
                // New column:
                ind++;
                finalvals[ind] = stitchedvals[reorderingvector[i]];
                finalcolindices[ind] = stitchedcolindices[reorderingvector[i]];
            }
        }
        else
        {
            // New row:
            row++; ind++;
        
            // If empty row move to the next not empty one:
            int newrow = stitchedrowindices[reorderingvector[i]];
            while (row < newrow) { finalrowindices[row] = ind; row++; }
            
            finalvals[ind] = stitchedvals[reorderingvector[i]];
            finalcolindices[ind] = stitchedcolindices[reorderingvector[i]];
            finalrowindices[row] = ind;
        }
    }
    while (row < countrows()) { row++; finalrowindices[row] = nnz; }
    
    delete[] stitchedrowindices;
    delete[] stitchedcolindices;
    delete[] stitchedvals;
    delete[] reorderingvector;
    
    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, countrows(), countcolumns(), finalrowindices, finalcolindices, finalvals, &mymat);
    
    MatAssemblyBegin(mymat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mymat, MAT_FINAL_ASSEMBLY);
}

void rawmat::removelastfragment(void)
{
    accumulatedrowindices.pop_back();
    accumulatedcolindices.pop_back();
    accumulatedvals.pop_back();
}

void rawmat::clearfragments(void)
{
    accumulatedrowindices = {};
    accumulatedcolindices = {};
    accumulatedvals = {};
}

void rawmat::print(void)
{            
    std::cout << std::endl;
    if (mydofmanager == NULL)
        std::cout << "Undefined matrix" << std::endl;
    else
    {
        if (countrows() > 0)
            MatView(mymat, PETSC_VIEWER_STDOUT_SELF);
        else
            std::cout << "Empty matrix" << std::endl;
    }
    std::cout << std::endl;
}

std::shared_ptr<rawmat> rawmat::extractaccumulated(void)
{
    std::shared_ptr<rawmat> output(new rawmat(mydofmanager));
    
    output->accumulatedrowindices = accumulatedrowindices;
    output->accumulatedcolindices = accumulatedcolindices;
    output->accumulatedvals = accumulatedvals;
    
    return output;
}



