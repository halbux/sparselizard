#include "disjointregions.h"
#include "universe.h"


int disjointregions::count(void)
{
    return elementtypenumbers.size();
}

int disjointregions::countelements(int disjointregionnumber)
{
    return rangeend[disjointregionnumber] - rangebegin[disjointregionnumber] + 1;
}

int disjointregions::add(int elementtypenumber, std::vector<bool>& physicalregionsincludingit)
{
    // Check if the disjoint region is already defined (loop is last to first for speedup):
    for (int disjreg = disjointregionsdefinition.size()-1; disjreg >= 0; disjreg--)
    {
        if (elementtypenumber == elementtypenumbers[disjreg] && disjointregionsdefinition[disjreg] == physicalregionsincludingit)
            return disjreg;
    }
    // In case it is not yet defined:
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

bool disjointregions::isinphysicalregion(int disjointregionnumber, int physicalregionindex)
{
    return disjointregionsdefinition[disjointregionnumber][physicalregionindex];
}

std::vector<int> disjointregions::get(int dim)
{
	int numdisjregs = 0;
	for (int i = 0; i < count(); i++)
	{
		if (getelementdimension(i) == dim)
			numdisjregs++;
	}
	std::vector<int> disjregs(numdisjregs);

	int index = 0;
	for (int i = 0; i < count(); i++)
	{
		if (getelementdimension(i) == dim)
		{
			disjregs[index] = i;
			index++;
		}
	}	
	return disjregs;
}

std::vector<int> disjointregions::getphysicalregions(int disjreg)
{
	std::vector<int> physregtranslation = universe::mymesh->getphysicalregions()->getallnumbers();

	int numphysregs = 0;
	for (int i = 0; i < disjointregionsdefinition[disjreg].size(); i++)
	{
		if (disjointregionsdefinition[disjreg][i])
			numphysregs++;
	}

	std::vector<int> physregs(numphysregs);

	int index = 0;
	for (int i = 0; i < disjointregionsdefinition[disjreg].size(); i++)
	{
		if (disjointregionsdefinition[disjreg][i])
		{
			physregs[index] = physregtranslation[i];
			index++;
		}
	}	
	return physregs;
}



