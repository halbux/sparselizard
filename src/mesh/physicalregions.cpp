#include "physicalregions.h"
#include "universe.h"


physicalregions::physicalregions(disjointregions& inputdisjointregions)
{
    mydisjointregions = &inputdisjointregions;
}

int physicalregions::createunion(std::vector<int> input, bool createifexisting)
{
    std::vector<int> disjregs = {};
    for (int i = 0; i < input.size(); i++)
    {
        // Get all disjoint regions with -1:
        std::vector<int> disjregsinthisphysreg = get(input[i])->getdisjointregions(-1);
        for (int j = 0; j < disjregsinthisphysreg.size(); j++)
            disjregs.push_back(disjregsinthisphysreg[j]);
    }

    if (createifexisting == false)
    {
        int existingnumber = find(disjregs);
        if (existingnumber >= 0)
            return existingnumber;
    }
    
    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(disjregs);
    
    return newphysregnum;
}

int physicalregions::createintersection(std::vector<int> input, bool createifexisting)
{
    std::vector<int> disjregs = {};
    for (int i = 0; i < input.size(); i++)
    {
        // Get all disjoint regions with -1:
        std::vector<int> disjregsinthisphysreg = get(input[i])->getdisjointregions(-1);
        
        if (i > 0)
            disjregs = gentools::intersect(disjregs, disjregsinthisphysreg);
        else
            disjregs = disjregsinthisphysreg;
    }
    
    if (createifexisting == false)
    {
        int existingnumber = find(disjregs);
        if (existingnumber >= 0)
            return existingnumber;
    }
    
    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(disjregs);
    
    return newphysregnum;
}

int physicalregions::createunionofall(bool createifexisting)
{
    std::vector<int> disjregs = gentools::getequallyspaced(0, 1, mydisjointregions->count());

    if (createifexisting == false)
    {
        int existingnumber = find(disjregs);
        if (existingnumber >= 0)
            return existingnumber;
    }
    
    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(disjregs);
    
    return newphysregnum;
}

int physicalregions::createfromdisjointregionlist(std::vector<int> drs)
{
    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(drs);
    
    return newphysregnum;
}

void physicalregions::definewithdisjointregions(void)
{
    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
        myphysicalregions[i]->definewithdisjointregions();
}

int physicalregions::getmaxphysicalregionnumber(void)
{
    return *std::max_element(myphysicalregionnumbers.begin(), myphysicalregionnumbers.end());
}

physicalregion* physicalregions::get(int physicalregionnumber)
{        
    // Try to find the physical region number in 'myphysicalregionnumbers':
    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
    {
        if (myphysicalregionnumbers[i] == physicalregionnumber)
            return myphysicalregions[i].get();
    }
    
    // If it could not be found append it:
    std::shared_ptr<physicalregion> newphysicalregion(new physicalregion(*mydisjointregions, *this, physicalregionnumber));
    myphysicalregions.push_back(newphysicalregion);
    myphysicalregionnumbers.push_back(physicalregionnumber);
    
    return (myphysicalregions[myphysicalregions.size()-1]).get();
}

physicalregion* physicalregions::getatindex(int physicalregionindex)
{
    return (myphysicalregions[physicalregionindex]).get();
}

int physicalregions::count(int dim)
{
    if (dim == -1)
        return myphysicalregionnumbers.size();
    else
    {
        int num = 0;
        for (int i = 0; i < myphysicalregionnumbers.size(); i++)
        {
            if (myphysicalregions[i]->getelementdimension() == dim)
                num++;
        }
        return num;
    }    
}

int physicalregions::countelements(void)
{
    int numelem = 0;
    for (int i = 0; i < myphysicalregions.size(); i++)
        numelem += myphysicalregions[i]->countelements();
    
    return numelem;
}

std::vector<int> physicalregions::getallnumbers(int dim)
{
    if (dim == -1)
        return myphysicalregionnumbers;
    else
    {
        std::vector<int> out(count(dim));
        int index = 0;
        for (int i = 0; i < myphysicalregions.size(); i++)
        {
            if (myphysicalregions[i]->getelementdimension() == dim)
            {   
                out[index] = myphysicalregionnumbers[i];
                index++;
            }
        }
        return out;
    }
}

int physicalregions::getnumber(int physicalregionindex)
{
    return myphysicalregionnumbers[physicalregionindex];
}

int physicalregions::getindex(int physicalregionnumber)
{
    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
    {
        if (myphysicalregionnumbers[i] == physicalregionnumber)
            return i;
    }
    return -1;
}

int physicalregions::find(std::vector<int>& disjregsinphysreg)
{
    std::vector<bool> argdef(mydisjointregions->count(), false);
    for (int i = 0; i < disjregsinphysreg.size(); i++)
        argdef[disjregsinphysreg[i]] = true;

    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
    {
        if (myphysicalregions[i]->getdefinition() == argdef)
            return myphysicalregionnumbers[i];
    }
    return -1;
}

void physicalregions::inphysicalregions(int elementtypenumber, int totalnumelemsintype, std::vector<int>& addresses, std::vector<int>& prs)
{
    element myelem(elementtypenumber);
    int dim = myelem.getelementdimension();

    std::vector<int> numprinelems(totalnumelemsintype, 0);
    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
    {
        if (myphysicalregions[i]->getelementdimension() != dim)
            continue;
        std::vector<int>* ellist = &(myphysicalregions[i]->getelementlist()->at(elementtypenumber));
        for (int j = 0; j < ellist->size(); j++)
            numprinelems[ellist->at(j)]++;
    }

    addresses = std::vector<int>(totalnumelemsintype+1, 0);
    for (int i = 1; i <= totalnumelemsintype; i++)
        addresses[i] = addresses[i-1] + numprinelems[i-1];

    // Populate:
    prs = std::vector<int>(addresses[totalnumelemsintype]);
    std::vector<int> ind(totalnumelemsintype, 0);
    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
    {
        if (myphysicalregions[i]->getelementdimension() != dim)
            continue;
        int curpr = myphysicalregionnumbers[i];
        std::vector<int>* ellist = &(myphysicalregions[i]->getelementlist()->at(elementtypenumber));
        for (int j = 0; j < ellist->size(); j++)
        {
            int curel = ellist->at(j);
            prs[addresses[curel]+ind[curel]] = curpr;
            ind[curel]++;
        }
    }
}

void physicalregions::remove(std::vector<int> toremove, bool ispartofdisjregstructure)
{
    // Tag regions to remove:
    std::vector<bool> istoremove(myphysicalregionnumbers.size(), false);

    for (int i = 0; i < toremove.size(); i++)
    {
        int curindex = getindex(toremove[i]);
        if (curindex != -1)
            istoremove[curindex] = true;
    }

    int index = 0;
    for (int i = 0; i < myphysicalregionnumbers.size(); i++)
    {
        if (not(istoremove[i]))
        {
            myphysicalregions[index] = myphysicalregions[i];
            myphysicalregionnumbers[index] = myphysicalregionnumbers[i];
            
            index++;
        }
    }
    myphysicalregions.resize(index); 
    myphysicalregionnumbers.resize(index);
    
    if (ispartofdisjregstructure)
        mydisjointregions->removephysicalregions(istoremove);
}

void physicalregions::errorundefined(std::vector<int> physregs)
{
    for (int i = 0; i < physregs.size(); i++)
    {
        if (getindex(physregs[i]) == -1)
        {
            std::cout << "Error in 'physicalregions' object: physical region number " << physregs[i] << " is not defined" << std::endl;
            abort();
        }
    }
}

void physicalregions::errornotsamedim(std::vector<int> physregs)
{
    if (physregs.size() <= 1)
        return;

    int dim = get(physregs[0])->getelementdimension();

    for (int i = 1; i < physregs.size(); i++)
    {
        int curdim = get(physregs[i])->getelementdimension();
        if (dim != curdim)
        {
            std::cout << "Error in 'physicalregions' object: expected physical regions of same dimension" << std::endl;
            abort();
        }
    }
}

void physicalregions::copy(disjointregions* drs, physicalregions* target)
{
    *target = *this;
    
    for (int i = 0; i < myphysicalregions.size(); i++)
        target->myphysicalregions[i] = myphysicalregions[i]->copy(target, drs);
    
    target->mydisjointregions = drs;
}

