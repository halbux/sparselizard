#include "disjointregions.h"


int disjointregions::count(void)
{
    return elementtypenumbers.size();
}

int disjointregions::countelements(int disjointregionnumber)
{
    return rangeend[disjointregionnumber] - rangebegin[disjointregionnumber] + 1;
}

int disjointregions::append(int elementtypenumber, std::vector<bool>& physicalregionsincludingit)
{
    disjointregionsdefinition.push_back(physicalregionsincludingit);
    elementtypenumbers.push_back(elementtypenumber);
    return disjointregionsdefinition.size() - 1;
}

void disjointregions::setrangebegin(int disjointregionnumber, int startrange)
{
    if (disjointregionnumber >= 0)
    {
        if (rangebegin.size() != elementtypenumbers.size())
            rangebegin.resize(elementtypenumbers.size());
        rangebegin[disjointregionnumber] = startrange;
    }
}

void disjointregions::setrangeend(int disjointregionnumber, int endrange)
{
    if (disjointregionnumber >= 0)
    {
        if (rangeend.size() != elementtypenumbers.size())
            rangeend.resize(elementtypenumbers.size());
        rangeend[disjointregionnumber] = endrange;
    }
}

int disjointregions::getrangebegin(int disjointregionnumber)
{
    return rangebegin[disjointregionnumber];
}

int disjointregions::getrangeend(int disjointregionnumber)
{
    return rangeend[disjointregionnumber];
}

int disjointregions::getelementtypenumber(int disjointregionnumber)
{
    return elementtypenumbers[disjointregionnumber];
}

int disjointregions::getelementdimension(int disjointregionnumber)
{
    element myelement(elementtypenumbers[disjointregionnumber]);
    return myelement.getelementdimension();
}

std::vector<int> disjointregions::getindim(int dim)
{
    std::vector<int> output(elementtypenumbers.size());
    int num = 0;
    for (int d = 0; d < output.size(); d++)
    {
        if (getelementdimension(d) == dim)
        {
            output[num] = d;
            num++;
        }
    }
    output.resize(num);
    
    return output;
}

std::vector<int> disjointregions::getintype(int elementtypenumber)
{
    std::vector<int> output(elementtypenumbers.size());
    int num = 0;
    for (int d = 0; d < output.size(); d++)
    {
        if (getelementtypenumber(d) == elementtypenumber)
        {
            output[num] = d;
            num++;
        }
    }
    output.resize(num);
    
    return output;
}

bool disjointregions::isinphysicalregion(int disjointregionnumber, int physicalregionindex)
{
    return disjointregionsdefinition[disjointregionnumber][physicalregionindex];
}

void disjointregions::removephysicalregions(std::vector<bool> istoremove)
{
    for (int i = 0; i < disjointregionsdefinition.size(); i++)
    {
        int index = 0;
        for (int j = 0; j < disjointregionsdefinition[i].size(); j++)
        {
            if (not(istoremove[j]))
            {
                disjointregionsdefinition[i][index] = disjointregionsdefinition[i][j];
                index++;    
            }
        }
        disjointregionsdefinition[i].resize(index);
    }
}

void disjointregions::clear(void)
{
    rangebegin = {};
    rangeend = {};
    disjointregionsdefinition = {};
    elementtypenumbers = {};
}

