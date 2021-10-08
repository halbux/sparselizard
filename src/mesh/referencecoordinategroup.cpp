#include "referencecoordinategroup.h"


referencecoordinategroup::referencecoordinategroup(std::vector<double>& coords)
{
    mycoordgroup = coordinategroup(coords);
}

referencecoordinategroup::referencecoordinategroup(std::vector<int>& elems, std::vector<double>& refcoords)
{
    myinputelems = elems;
    myinputrefcoords = refcoords;
}

void referencecoordinategroup::evalat(std::vector<int> inputdisjregs)
{
    int numcoords = mycoordgroup.countcoordinates();
    
    std::vector<int> elems(numcoords,-1);
    std::vector<int> coordnums(numcoords);
    std::iota(coordnums.begin(), coordnums.end(), 0);
    std::vector<double> kietaphis(3*numcoords,0.0);
    
    for (int d = 0; d < inputdisjregs.size(); d++)
        myalgorithm::getreferencecoordinates(mycoordgroup, inputdisjregs[d], elems, kietaphis);
        
    evalat(elems, kietaphis, coordnums);  
}

void referencecoordinategroup::evalat(int elemtypenum)
{
    int numintype = 0;
    for (int i = 0; i < myinputelems.size()/2; i++)
    {
        if (myinputelems[2*i+0] == elemtypenum)
            numintype++;
    }
    
    std::vector<int> elems(numintype);
    std::vector<double> kietaphis(3*numintype);
    std::vector<int> coordnums(numintype);

    int index = 0;
    for (int i = 0; i < myinputelems.size()/2; i++)
    {
        if (myinputelems[2*i+0] == elemtypenum)
        {
            elems[index] = myinputelems[2*i+1];
            kietaphis[3*index+0] = myinputrefcoords[3*i+0];
            kietaphis[3*index+1] = myinputrefcoords[3*i+1];
            kietaphis[3*index+2] = myinputrefcoords[3*i+2];
            coordnums[index] = i;
            
            index++;
        }
    }
    
    evalat(elems, kietaphis, coordnums);  
}
    
void referencecoordinategroup::evalat(std::vector<int>& elems, std::vector<double>& kietaphis, std::vector<int>& coordnums)
{
    int numcoords = elems.size();
    
    // Remove the negative elements (not found):
    int numpositive = 0;
    for (int i = 0; i < numcoords; i++)
    {
        if (elems[i] >= 0)
            numpositive++;
    }
    if (numpositive == 0)
        return;

    int curindex = 0;
    for (int i = 0; i < numcoords; i++)
    {   
        if (elems[i] >= 0)
        {
            elems[curindex] = elems[i];
            coordnums[curindex] = coordnums[i];
            kietaphis[3*curindex+0] = kietaphis[3*i+0];
            kietaphis[3*curindex+1] = kietaphis[3*i+1];
            kietaphis[3*curindex+2] = kietaphis[3*i+2];
            curindex++;
        }
    }
    elems.resize(numpositive);
    coordnums.resize(numpositive);
    kietaphis.resize(3*numpositive);
    
    // Sort according to the element numbers:
    std::vector<int> reorderingvector;
    myalgorithm::stablesort(elems, reorderingvector);
    std::vector<int> elemscp = elems;
    std::vector<int> coordnumscp = coordnums;
    std::vector<double> kietaphiscp = kietaphis;
    for (int i = 0; i < elems.size(); i++)
    {
        elems[i] = elemscp[reorderingvector[i]];
        coordnums[i] = coordnumscp[reorderingvector[i]];
        kietaphis[3*i+0] = kietaphiscp[3*reorderingvector[i]+0];
        kietaphis[3*i+1] = kietaphiscp[3*reorderingvector[i]+1];
        kietaphis[3*i+2] = kietaphiscp[3*reorderingvector[i]+2];
    }

    // Count the number of elements with a given number of ref coords:
    int maxnumrefcoords = 1; int curnumrefcoords = 1;
    for (int i = 1; i < elems.size(); i++)
    {
        if (elems[i] == elems[i-1])
        {
            curnumrefcoords++;
            if (maxnumrefcoords < curnumrefcoords)
                maxnumrefcoords = curnumrefcoords;
        }
        else
            curnumrefcoords = 1;
    }
    std::vector<int> numelemswithnumrefcoords(maxnumrefcoords+1,0);
    numelemswithnumrefcoords[1] = 1;
    curnumrefcoords = 1;
    for (int i = 1; i < elems.size(); i++)
    {
        if (elems[i] == elems[i-1])
        {
            numelemswithnumrefcoords[curnumrefcoords]--;
            numelemswithnumrefcoords[curnumrefcoords+1]++;
            curnumrefcoords++;
        }
        else
        {
            numelemswithnumrefcoords[1]++;
            curnumrefcoords = 1;
        }
    }
    
    // Split the elements into number of reference coordinates:
    myelems.resize(maxnumrefcoords+1);
    mycoordnums.resize(maxnumrefcoords+1);
    mykietaphis.resize(maxnumrefcoords+1);
    for (int n = 1; n <= maxnumrefcoords; n++)
    {
        int numelems = numelemswithnumrefcoords[n];
        myelems[n].resize(numelems*n);
        mycoordnums[n].resize(numelems*n);
        mykietaphis[n].resize(3*numelems*n);
    }
    std::vector<int> curindexes(maxnumrefcoords+1,0);
    int first = 0;
    for (int i = 1; i <= elems.size(); i++)
    {
        if (i < elems.size() && elems[i] == elems[i-1])
            continue;
        else
        {
            // Number of ref coords for the current element:
            int n = i-first;
    
            for (int j = 0; j < n; j++)
            {
                myelems[n][curindexes[n]] = elems[first+j];
                mycoordnums[n][curindexes[n]] = coordnums[first+j];
                mykietaphis[n][3*curindexes[n]+0] = kietaphis[3*(first+j)+0];
                mykietaphis[n][3*curindexes[n]+1] = kietaphis[3*(first+j)+1];
                mykietaphis[n][3*curindexes[n]+2] = kietaphis[3*(first+j)+2];
                
                curindexes[n]++;
            }
            first += n;
        }
    }
    
    // Sort the containers according to the reference coordinates:
    for (int n = 1; n <= maxnumrefcoords; n++)
    {
        // Sort the reference coordinates in each element:
        std::vector<int> reordvec;
        myalgorithm::stablecoordinatesort({noisethreshold,noisethreshold,noisethreshold}, myelems[n], mykietaphis[n], reordvec);

        std::vector<int> myelemsreordered(myelems[n].size());
        std::vector<int> mycoordnumsreordered(myelems[n].size());
        std::vector<double> mykietaphisreordered(3*myelems[n].size());
        for (int i = 0; i < reordvec.size(); i++)
        {
            myelemsreordered[i] = myelems[n][reordvec[i]];
            mycoordnumsreordered[i] = mycoordnums[n][reordvec[i]];
            mykietaphisreordered[3*i+0] = mykietaphis[n][3*reordvec[i]+0];
            mykietaphisreordered[3*i+1] = mykietaphis[n][3*reordvec[i]+1];
            mykietaphisreordered[3*i+2] = mykietaphis[n][3*reordvec[i]+2];
        }
        // Sort by blocks of n coordinates:
        myalgorithm::stablesort(noisethreshold, mykietaphisreordered, reordvec, 3*n);
    
        for (int i = 0; i < reordvec.size(); i++)
        {
            for (int j = 0; j < n; j++)
            {
                myelems[n][n*i+j] = myelemsreordered[n*reordvec[i]+j];
                mycoordnums[n][n*i+j] = mycoordnumsreordered[n*reordvec[i]+j];
                mykietaphis[n][3*(n*i+j)+0] = mykietaphisreordered[3*(n*reordvec[i]+j)+0];
                mykietaphis[n][3*(n*i+j)+1] = mykietaphisreordered[3*(n*reordvec[i]+j)+1];
                mykietaphis[n][3*(n*i+j)+2] = mykietaphisreordered[3*(n*reordvec[i]+j)+2];
            }
        }
    }
    
    myrangebegin = 0; mynumrefcoords = 1;
}

bool referencecoordinategroup::next(void)
{
    // Skip empty blocks:
    while (mynumrefcoords < myelems.size() && myelems[mynumrefcoords].size() == 0)
        mynumrefcoords++;
        
    if (mynumrefcoords >= myelems.size())
        return false;
        
    std::vector<int>* curels = &(myelems[mynumrefcoords]);
    std::vector<int>* curcoords = &(mycoordnums[mynumrefcoords]);
    std::vector<double>* curkep = &(mykietaphis[mynumrefcoords]);
        
    // Get all reference coordinates in the first element:
    mycurrefcoords.resize(3*mynumrefcoords);
    for (int i = 0; i < 3*mynumrefcoords; i++)
        mycurrefcoords[i] = curkep->at(3*myrangebegin+i);
        
    // Find all consecutive elements with close enough reference coordinates:
    int curelemindex = myrangebegin+mynumrefcoords; bool continueit = true;
    while (curelemindex < curels->size() && continueit)
    {
        for (int i = 0; i < mynumrefcoords; i++)
        {
            double kiref = mycurrefcoords[3*i+0]; double etaref = mycurrefcoords[3*i+1]; double phiref = mycurrefcoords[3*i+2];
            double kicur = curkep->at(3*(curelemindex+i)+0); double etacur = curkep->at(3*(curelemindex+i)+1); double phicur = curkep->at(3*(curelemindex+i)+2);
            if (std::abs(kiref-kicur) > noisethreshold || std::abs(etaref-etacur) > noisethreshold || std::abs(phiref-phicur) > noisethreshold)
            {
                curelemindex -= mynumrefcoords;
                continueit = false;
                break;
            }
        }
        curelemindex += mynumrefcoords;
    }
    
    // Get the list of elements of same reference coordinates:
    int numelemsinrange = (curelemindex-myrangebegin)/mynumrefcoords;
    mycurelems.resize(numelemsinrange);
    for (int i = 0; i < numelemsinrange; i++)
        mycurelems[i] = curels->at(myrangebegin+mynumrefcoords*i);
    
    // Get the numbers of the reference coordinates:
    mycurcoordnums.resize(numelemsinrange*mynumrefcoords);
    for (int i = 0; i < mycurcoordnums.size(); i++)
        mycurcoordnums[i] = curcoords->at(myrangebegin+i);
    


    myrangebegin = curelemindex;
    
    if (myrangebegin == curels->size())
    {
        mynumrefcoords++;
        myrangebegin = 0;
    }
    
    
    return true;
}

