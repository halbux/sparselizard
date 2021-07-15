#include "rawmat.h"
#include <thread>


rawmat::rawmat(std::shared_ptr<dofmanager> dofmngr)
{
    mydofmanager = dofmngr;
    
    mymeshnumber = universe::mymesh->getmeshnumber();
}

rawmat::rawmat(std::shared_ptr<dofmanager> dofmngr, Mat inA, Mat inD, intdensematrix inAinds, intdensematrix inDinds)
{
    mydofmanager = dofmngr;
    
    Amat = inA;
    Dmat = inD;
    Ainds = inAinds;
    Dinds = inDinds;
    
    mymeshnumber = universe::mymesh->getmeshnumber();
}
        
rawmat::~rawmat(void) 
{
    // Avoid crashes when destroy is called after PetscFinalize (not allowed).
    PetscBool ispetscinitialized;
    PetscInitialized(&ispetscinitialized);

    if (ispetscinitialized == PETSC_TRUE)
    {
        if (Amat != PETSC_NULL)
            MatDestroy(&Amat);
        if (Dmat != PETSC_NULL)
            MatDestroy(&Dmat);
        if (isitfactored) 
            KSPDestroy(&myksp);
    }
}

long long int rawmat::countrows(void) 
{ 
    if (mydofmanager == NULL)
        return 0;
    else
        return mydofmanager->countdofs();
}

long long int rawmat::countcolumns(void) 
{ 
    if (mydofmanager == NULL)
        return 0;
    else
        return mydofmanager->countdofs();
}

void rawmat::accumulate(intdensematrix rowadresses, intdensematrix coladresses, densematrix vals)
{
    if (vals.count() > 0)
    {
        accumulatedrowindices.push_back(rowadresses);
        accumulatedcolindices.push_back(coladresses);
        accumulatedvals.push_back(vals);
    }
}

void processrows(int firstrow, int lastrow, int* maxnnzinrows, long long int* adsofrows, std::pair<int, double>* valsptr, std::vector<bool>* isconstrained, int* nnzApart, int* nnzDpart)
{
    // Avoid cache line invalidation:
    int curnnzA = 0, curnnzD = 0;
    
    for (int i = firstrow; i <= lastrow; i++)
    {
        int maxnnzinrow = maxnnzinrows[i];
        std::pair<int, double>* curvalsptr = valsptr + adsofrows[i];
        
        std::vector<int> curcols(maxnnzinrow);
        std::vector<double> curvals(maxnnzinrow);
        for (int j = 0; j < maxnnzinrow; j++)
        {
            curcols[j] = curvalsptr[j].first;
            curvals[j] = curvalsptr[j].second;
        }
        
        // Single threaded sort:
        std::vector<int> reorderingvector(maxnnzinrow);
        std::iota(reorderingvector.begin(), reorderingvector.end(), 0);
        std::sort(reorderingvector.begin(), reorderingvector.end(), [&](int elem1, int elem2)
        { 
            if (curcols[elem1] < curcols[elem2])
                return true;
            if (curcols[elem1] > curcols[elem2])
                return false;
            return elem1 < elem2;
        });
    
        int pc = -1;
        int ind = -1;
        for (int j = 0; j < maxnnzinrow; j++)
        {
            int cc = curcols[reorderingvector[j]];
            double cv = curvals[reorderingvector[j]];
            
            if (cc == pc)
                curvalsptr[ind].second += cv;
            else
            {
                ind++;
                curvalsptr[ind].first = cc;
                curvalsptr[ind].second = cv;

                if (isconstrained->at(cc))
                    curnnzD++;
                else
                    curnnzA++;
            }
                
            pc = cc;
        }
        maxnnzinrows[i] = ind+1;
    }
    
    *nnzApart = curnnzA;
    *nnzDpart = curnnzD;
}

void rawmat::process(std::vector<bool>& isconstrained)
{
    int ndofs = countrows();
    
    // Create Ainds and Dinds:
    std::vector<int> renumtolocalindex;
    myalgorithm::findtruefalse(isconstrained, Dinds, Ainds, renumtolocalindex);
    
    // Get an upper bound on the number of nonzeros in each row:
    long long int maxnnz = 0;
    std::vector<int> maxnnzinrows(ndofs, 0);
    for (int i = 0; i < accumulatedvals.size(); i++)
    {
        int* accumulatedrowindicesptr = accumulatedrowindices[i].getvalues();
        int* accumulatedcolindicesptr = accumulatedcolindices[i].getvalues();
        
        int nr = accumulatedvals[i].countrows();
        int nc = accumulatedvals[i].countcolumns();
        int ntr = accumulatedrowindices[i].countrows();
        int ndr = accumulatedcolindices[i].countrows();
        
        for (int r = 0; r < nr; r++)
        {
            int ctr = r, cdr = r;
            if (ntr != nr || ndr != nr)
            {
                ctr = r/ndr;
                cdr = r%ndr;
            }
        
            for (long long int c = 0; c < nc; c++)
            {
                int cr = accumulatedrowindicesptr[ctr*nc+c];
                int cc = accumulatedcolindicesptr[cdr*nc+c];
                
                if (cr >= 0 && cc >= 0 && isconstrained[cr] == false)
                {
                    maxnnzinrows[cr]++;
                    maxnnz++;
                }
            }
        } 
    }

    // Create a vector for direct addressing:
    std::vector<long long int> adsofrows(ndofs, 0);
    for (int i = 1; i < ndofs; i++)
        adsofrows[i] = adsofrows[i-1] + maxnnzinrows[i-1];
        
    // Collect the values for each row:
    std::vector<std::pair<int, double>> valspairs(maxnnz);
    std::pair<int, double>* valsptr = valspairs.data();
    
    std::vector<int> indexinrow(ndofs, 0);
    for (int i = 0; i < accumulatedvals.size(); i++)
    {
        int* accumulatedrowindicesptr = accumulatedrowindices[i].getvalues();
        int* accumulatedcolindicesptr = accumulatedcolindices[i].getvalues();
        double* accumulatedvalsptr = accumulatedvals[i].getvalues();
        
        int nr = accumulatedvals[i].countrows();
        int nc = accumulatedvals[i].countcolumns();
        int ntr = accumulatedrowindices[i].countrows();
        int ndr = accumulatedcolindices[i].countrows();
        
        for (int r = 0; r < nr; r++)
        {
            int ctr = r, cdr = r;
            if (ntr != nr || ndr != nr)
            {
                ctr = r/ndr;
                cdr = r%ndr;
            }
        
            for (long long int c = 0; c < nc; c++)
            {
                int cr = accumulatedrowindicesptr[ctr*nc+c];
                int cc = accumulatedcolindicesptr[cdr*nc+c];
                
                double cv = accumulatedvalsptr[r*nc+c];
                    
                if (cr >= 0 && cc >= 0 && isconstrained[cr] == false)
                {
                    long long int curind = adsofrows[cr] + indexinrow[cr];
                    valsptr[curind].first = cc;
                    valsptr[curind].second = cv;
                    
                    indexinrow[cr]++;
                }
            }
        } 
    }

    // Multithreaded sort and duplicate removal:
    int numthreadstouse = std::min(ndofs/10000+1, universe::getmaxnumthreads()); // require a min num dofs per thread

    std::vector<int> nnzAparts(numthreadstouse, 0), nnzDparts(numthreadstouse, 0);    
    if (numthreadstouse > 1)
    {
        std::vector<std::thread> threadobjs(numthreadstouse);
        int rowchunksize = ndofs/numthreadstouse+1;

        for (int t = 0; t < numthreadstouse; t++)
            threadobjs[t] = std::thread(processrows, t*rowchunksize, std::min((t+1)*rowchunksize-1, ndofs-1), maxnnzinrows.data(), adsofrows.data(), valsptr, &isconstrained, &nnzAparts[t], &nnzDparts[t]);
        
        for (int t = 0; t < numthreadstouse; t++)
            threadobjs[t].join();
    }
    else
        processrows(0, ndofs-1, maxnnzinrows.data(), adsofrows.data(), valsptr, &isconstrained, &nnzAparts[0], &nnzDparts[0]);

    nnzA = myalgorithm::sum(nnzAparts);
    nnzD = myalgorithm::sum(nnzDparts);

    // Create A and D:
    Arows = intdensematrix(Ainds.count()+1,1);
    Acols = intdensematrix(nnzA, 1);
    Avals = densematrix(nnzA, 1);
    int* Arowsptr = Arows.getvalues();
    int* Acolsptr = Acols.getvalues();
    double* Avalsptr = Avals.getvalues();
    
    Drows = intdensematrix(Ainds.count()+1,1);
    Dcols = intdensematrix(nnzD, 1);
    Dvals = densematrix(nnzD, 1);
    int* Drowsptr = Drows.getvalues();
    int* Dcolsptr = Dcols.getvalues();
    double* Dvalsptr = Dvals.getvalues();
    
    Arowsptr[0] = 0;
    Drowsptr[0] = 0;

    int index = 0;
    for (int i = 0; i < ndofs; i++)
    {
        if (isconstrained[i])
            continue;
    
        Arowsptr[index+1] = Arowsptr[index];
        Drowsptr[index+1] = Drowsptr[index];
        
        std::pair<int, double>* curvalsptr = valsptr + adsofrows[i];        

        for (int j = 0; j < maxnnzinrows[i]; j++)
        {
            int cc = curvalsptr[j].first;
            double cv = curvalsptr[j].second;
            
            if (isconstrained[cc])
            {
                Dcolsptr[Drowsptr[index+1]] = renumtolocalindex[cc];
                Dvalsptr[Drowsptr[index+1]] = cv;
                Drowsptr[index+1]++;
            }
            else
            {
                Acolsptr[Arowsptr[index+1]] = renumtolocalindex[cc];
                Avalsptr[Arowsptr[index+1]] = cv;
                Arowsptr[index+1]++;
            }
        }
        index++;
    }

    // Create the petsc objects:
    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, Ainds.count(), Ainds.count(), Arowsptr, Acolsptr, Avalsptr, &Amat);
    MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY);

    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, Ainds.count(), Dinds.count(), Drowsptr, Dcolsptr, Dvalsptr, &Dmat);
    MatAssemblyBegin(Dmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Dmat, MAT_FINAL_ASSEMBLY);
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
        {
            std::cout << "A block " << Ainds.count() << "x" << Ainds.count() << ":" << std::endl << std::endl;
            MatView(Amat, PETSC_VIEWER_STDOUT_SELF);
            std::cout << std::endl << "D block " << Ainds.count() << "x" << Dinds.count() << ":" << std::endl << std::endl;
            MatView(Dmat, PETSC_VIEWER_STDOUT_SELF);
        }
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

intdensematrix rawmat::getainds(void)
{
    return Ainds;
}

intdensematrix rawmat::getdinds(void)
{
    return Dinds;
}

Mat rawmat::getapetsc(void)
{
    return Amat;
}

Mat rawmat::getdpetsc(void)
{
    return Dmat;
}

std::shared_ptr<dofmanager> rawmat::getdofmanager(void)
{
    return mydofmanager;
}

KSP* rawmat::getksp(void)
{
    return &myksp;
}

