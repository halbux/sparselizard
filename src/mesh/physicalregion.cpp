#include "physicalregion.h"


physicalregion::physicalregion(disjointregions& inputdisjointregions, physicalregions& inputphysicalregions, int physicalregionnumber)
{
    mydisjointregions = &inputdisjointregions;
    myphysicalregions = &inputphysicalregions;
    myphysicalregionnumber = physicalregionnumber;
}

int physicalregion::getnumber(void)
{
    return myphysicalregionnumber;
}

void physicalregion::addelement(int elementtypenumber, int elementnumber)
{
    element myelement(elementtypenumber);
    // Set once and for all the element dimension if not already done:
    if (myelementdimension == -1)
        myelementdimension = myelement.getelementdimension();
    
    if (myelementdimension != myelement.getelementdimension())
    {
        std::cout << "Error in 'physicalregion' object: trying to add a " << myelement.getelementdimension() << "D element to physical region " << myphysicalregionnumber << " which has " << myelementdimension << "D elements." << std::endl;
        std::cout << "There can only be a single element dimension per physical region." << std::endl;
        abort();
    }
    
    elementlist[elementtypenumber].push_back(elementnumber);
}

int physicalregion::countelements(void)
{
    std::vector<int> alldisjointregions = getdisjointregions();
    int numelem = 0;
    for (int i = 0; i < alldisjointregions.size(); i++)
        numelem += mydisjointregions->countelements(alldisjointregions[i]);
    return numelem;
}

int physicalregion::getelementdimension(void)
{
    return myelementdimension;
}

void physicalregion::definewithdisjointregions(void)
{   
    includesdisjointregion.resize(mydisjointregions->count());
    std::fill(includesdisjointregion.begin(), includesdisjointregion.end(), false);

    int prindex = myphysicalregions->getindex(myphysicalregionnumber);
    
    for (int i = 0; i < mydisjointregions->count(); i++)
        includesdisjointregion[i] = mydisjointregions->isinphysicalregion(i, prindex);
}

void physicalregion::setdisjointregions(std::vector<int> disjointregionlist)
{
    if (disjointregionlist.size() == 0)
    {
        std::cout << "Error in 'physicalregion' object: physical region cannot be empty" << std::endl;
        abort();
    }
    
    includesdisjointregion.resize(mydisjointregions->count());
    std::fill(includesdisjointregion.begin(), includesdisjointregion.end(), false);

    myelementdimension = mydisjointregions->getelementdimension(disjointregionlist[0]);
    for (int i = 0; i < disjointregionlist.size(); i++)
    {
        includesdisjointregion[disjointregionlist[i]] = true;
        if (myelementdimension < mydisjointregions->getelementdimension(disjointregionlist[i]))
            myelementdimension = mydisjointregions->getelementdimension(disjointregionlist[i]);
    }
}

std::vector<int> physicalregion::getdisjointregions(void)
{
    std::vector<int> disjointregionlist;
    for (int i = 0; i < includesdisjointregion.size(); i++)
    {
        // Only disjoint regions including elements of the same 
        // dimension as in the physical region are selected.
        if (includesdisjointregion[i] && myelementdimension == mydisjointregions->getelementdimension(i))
            disjointregionlist.push_back(i);
    }
    return disjointregionlist;
}

std::vector<int> physicalregion::getdisjointregions(int dim)
{
    std::vector<int> disjointregionlist;
    for (int i = 0; i < includesdisjointregion.size(); i++)
    {
        if (includesdisjointregion[i])
        {
            // If dim is -1 we take them all:
            if (mydisjointregions->getelementdimension(i) == dim || dim == -1)
                disjointregionlist.push_back(i);
        }
    }
    return disjointregionlist;
}

void physicalregion::renumberelements(int elementtypenumber, std::vector<int>& elementrenumbering)
{
    // Update 'elementlist'.
    for (int i = 0; i < elementlist[elementtypenumber].size(); i++)
        elementlist[elementtypenumber][i] = elementrenumbering[elementlist[elementtypenumber][i]];
}

void physicalregion::removeduplicatedelements(void)
{
    // Iterate on all element types:
    for (int i = 0; i <= 7; i++)
    {
        // Skip empty element types:
        if (elementlist[i].size() == 0)
            continue;

        // Sort every vector:
        sort(elementlist[i].begin(),elementlist[i].end());
        // Unique every vector and resize it to fit:
        std::vector<int>::iterator iter;
        iter = std::unique(elementlist[i].begin(),elementlist[i].end());
        elementlist[i].resize(std::distance(elementlist[i].begin(),iter));
    }
}

std::vector<std::vector<int>>* physicalregion::getelementlist(void)
{
    return &elementlist;
}


