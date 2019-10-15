#include "referencecoordinategroup.h"
#include "myalgorithm.h"


referencecoordinategroup::referencecoordinategroup(std::vector<double>& coords)
{
    mynumcoords = coords.size()/3;
    mycoordgroup = coordinategroup(coords);
}

void referencecoordinategroup::evalat(std::vector<int> inputdisjregs)
{
    mydisjregs = inputdisjregs;
    
    myelems = std::vector<int>(mynumcoords,-1);
    mykietaphis = std::vector<double>(3*mynumcoords,0.0);
    
    for (int d = 0; d < mydisjregs.size(); d++)
        myalgorithm::getreferencecoordinates(mycoordgroup, mydisjregs[d], myelems, mykietaphis);
        
    int myrangebegin = 0;
}

void referencecoordinategroup::group(void)
{
    // Sort first according to the element number then according to the reference coordinates:
    std::vector<int> myreorderingvector;
    myalgorithm::stablecoordinatesort({noisethreshold,noisethreshold,noisethreshold}, myelems, mykietaphis, myreorderingvector);
    
    std::vector<int> myelemscopy = myelems;
    std::vector<double> mykietaphiscopy = mykietaphis;
    for (int i = 0; i < mynumcoords; i++)
    {
        myelems[i] = myelemscopy[myreorderingvector[i]];
        mykietaphis[3*i+0] = mykietaphiscopy[3*myreorderingvector[i]+0];
        mykietaphis[3*i+1] = mykietaphiscopy[3*myreorderingvector[i]+1];
        mykietaphis[3*i+2] = mykietaphiscopy[3*myreorderingvector[i]+2];
    }
    myunorderingvector.resize(mynumcoords);
    for (int i = 0; i < mynumcoords; i++)
        myunorderingvector[myreorderingvector[i]] = i;
    
    // Move past the negative element numbers:
    while (myelems[myrangebegin] < 0)
        myrangebegin++;
}

bool referencecoordinategroup::next(void)
{
    if (myrangebegin >= mynumcoords)
        return false;
        
    // Get all reference coordinates in the first element:
    int index = myrangebegin+1;
    while (index < mynumcoords && myelems[index] == myelems[index-1])
        index++;
    int numrefcoords = index-myrangebegin;
    mycurrefcoords.resize(3*numrefcoords);
    for (int i = 0; i < 3*numrefcoords; i++)
        mycurrefcoords[i] = mykietaphis[myrangebegin+i];
    int curelem = myelems[myrangebegin];
        
    // Find all consecutive elements with close enough reference coordinates:
    int curelemindex = myrangebegin+numrefcoords; bool continueit = true;
    while (curelemindex < mynumcoords-(numrefcoords-1) && continueit)
    {
        for (int i = 0; i < numrefcoords; i++)
        {
            double kiref = mycurrefcoords[3*i+0]; double etaref = mycurrefcoords[3*i+1]; double phiref = mycurrefcoords[3*i+2];
            double kicur = mykietaphis[3*(curelemindex+i)+0]; double etacur = mykietaphis[3*(curelemindex+i)+1]; double phicur = mykietaphis[3*(curelemindex+i)+2];
            if (myelems[curelemindex+i] != curelem || std::abs(kiref-kicur) > noisethreshold || std::abs(etaref-etacur) > noisethreshold || std::abs(phiref-phicur) > noisethreshold)
            {
                curelemindex -= numrefcoords;
                continueit = false;
                break;
            }
        }
        curelemindex += numrefcoords;
    }
    
    // Get the list of elements of same reference coordinates:
    int numelemsinrange = (curelemindex-myrangebegin)/numrefcoords;
    mycurelems.resize(numelemsinrange);
    for (int i = 0; i < numelemsinrange; i++)
        mycurelems[i] = myelems[myrangebegin+numrefcoords*i];
    
    // Get the numbers of the reference coordinates:
    mycurcoordnums.resize(numelemsinrange*numrefcoords);
    for (int i = 0; i < mycurcoordnums.size(); i++)
        mycurcoordnums[i] = myunorderingvector[myrangebegin+i];
    

    myrangebegin = curelemindex;
    
    return true;
}

