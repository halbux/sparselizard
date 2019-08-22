#include "physicalregions.h"


physicalregions::physicalregions(disjointregions& inputdisjointregions)
{
    mydisjointregions = &inputdisjointregions;
}

int physicalregions::createunion(const std::vector<int> input)
{
    errorundefined(input);

    std::vector<int> disjregs = {};
    for (int i = 0; i < input.size(); i++)
    {
        // Get all disjoint regions with -1:
        std::vector<int> disjregsinthisphysreg = get(input[i])->getdisjointregions(-1);
        for (int j = 0; j < disjregsinthisphysreg.size(); j++)
            disjregs.push_back(disjregsinthisphysreg[j]);
    }
    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(disjregs);
    
    return newphysregnum;
}

int physicalregions::createintersection(const std::vector<int> input)
{
    errorundefined(input);

    std::vector<int> disjregs = {};
    for (int i = 0; i < input.size(); i++)
    {
        // Get all disjoint regions with -1:
        std::vector<int> disjregsinthisphysreg = get(input[i])->getdisjointregions(-1);
        
        if (i > 0)
            disjregs = myalgorithm::intersect(disjregs, disjregsinthisphysreg);
        else
            disjregs = disjregsinthisphysreg;
    }
    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(disjregs);
    
    return newphysregnum;
}

int physicalregions::createexclusion(int input, int toexclude)
{
    errorundefined({input,toexclude});

    // Get all input disjoint regions:
    std::vector<int> inputdisjregs = get(input)->getdisjointregions(-1);

    // Get all disjoint regions to exclude:
    std::vector<int> toexcludedisjregs = get(toexclude)->getdisjointregions(-1);

    // Remove the disjoint regions to exclude from the list of input disjoint regions:
    std::vector<int> disjregs = {};
    for (int i = 0; i < inputdisjregs.size(); i++)
    {
        bool toinclude = true;
        for (int j = 0; j < toexcludedisjregs.size(); j++)
        {
            if (inputdisjregs[i] == toexcludedisjregs[j])
            {
                toinclude = false;
                break;            
            }
        }
        if (toinclude)
            disjregs.push_back(inputdisjregs[i]);
    }

    int newphysregnum = getmaxphysicalregionnumber() + 1;
    
    physicalregion* newphysreg = get(newphysregnum);
    newphysreg->setdisjointregions(disjregs);
    
    return newphysregnum;
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
            return &myphysicalregions[i];
    }
    
    // If it could not be found append it:
    physicalregion newphysicalregion(*mydisjointregions, physicalregionnumber, myphysicalregionnumbers.size());
    myphysicalregions.push_back(newphysicalregion);
    myphysicalregionnumbers.push_back(physicalregionnumber);
    
    return &(myphysicalregions[myphysicalregions.size()-1]);
}

physicalregion* physicalregions::getatindex(int physicalregionindex)
{
    return &(myphysicalregions[physicalregionindex]);
}

int physicalregions::count(void)
{
    return myphysicalregionnumbers.size();
}

int physicalregions::countelements(void)
{
    int numelem = 0;
    for (int i = 0; i < myphysicalregions.size(); i++)
        numelem += myphysicalregions[i].countelements();
    
    return numelem;
}

std::vector<int> physicalregions::getallnumbers(void)
{
    return myphysicalregionnumbers;
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

