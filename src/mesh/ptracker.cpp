#include "ptracker.h"


ptracker::ptracker(std::vector<int> numelemspertype)
{
    elementrenumbering.resize(8);
    for (int i = 0; i < 8; i++)
        elementrenumbering[i] = myalgorithm::getequallyspaced(0, 1, numelemspertype[i]);
}

void ptracker::updatedisjointregions(disjointregions* input)
{
    mydisjointregions = *input;
}

void ptracker::updaterenumbering(std::vector<std::vector<int>>& renumber)
{
    for (int i = 0; i < 8; i++)
    {
        std::vector<int> oldrenum = elementrenumbering[i]; 
        for (int j = 0; j < elementrenumbering[i].size(); j++)
            elementrenumbering[i][j] = renumber[i][oldrenum[j]];
    }
}

void ptracker::getrenumbering(std::shared_ptr<ptracker> mt, std::vector<std::vector<int>>& renumbering)
{
    renumbering.resize(8);

    for (int i = 0; i < 8; i++)
    {
        renumbering[i].resize(elementrenumbering[i].size());
        if (mt != NULL)
        {
            for (int j = 0; j < elementrenumbering[i].size(); j++)
                renumbering[i][elementrenumbering[i][j]] = mt->elementrenumbering[i][j];
        }
        else
        {
            for (int j = 0; j < elementrenumbering[i].size(); j++)
                renumbering[i][elementrenumbering[i][j]] = j;
        }
    }
}

disjointregions* ptracker::getdisjointregions(void)
{
    return &mydisjointregions;
}

void ptracker::getindisjointregions(std::vector<std::vector<int>>& indisjregs)
{
    // Allocate the vectors:
    indisjregs.resize(8);
    for (int i = 0; i < 8; i++)
        indisjregs[i] = std::vector<int>(elementrenumbering[i].size(), -1); // Corner nodes will have -1
    
    int numdisjregs = mydisjointregions.count();
    for (int d = 0; d < numdisjregs; d++)
    {
        int typenum = mydisjointregions.getelementtypenumber(d);
        int rangebegin = mydisjointregions.getrangebegin(d);
        int numelems = mydisjointregions.countelements(d);
        
        for (int i = 0; i < numelems; i++)
            indisjregs[typenum][rangebegin+i] = d;
    }
}

void ptracker::print(void)
{
    for (int i = 0; i < elementrenumbering.size(); i++)
    {
        std::cout << std::endl << "Element type number " << i << ":" << std::endl << std::endl;
        for (int j = 0; j < elementrenumbering[i].size(); j++)
            std::cout << j << " --> " << elementrenumbering[i][j] << std::endl;
    }
    std::cout << std::endl;
}

