#include "meshtracker.h"


meshtracker::meshtracker(void) {}

void meshtracker::updatedisjointregions(disjointregions* input)
{
    mydisjointregions = *input;
}

void meshtracker::updaterenumbering(std::vector<std::vector<int>>& renumber)
{
    if (elementrenumbering.size() == 0)
        elementrenumbering = renumber;
    else
    {
        for (int i = 0; i < 8; i++)
        {
            std::vector<int> oldrenum = elementrenumbering[i]; 
            for (int j = 0; j < elementrenumbering[i].size(); j++)
                elementrenumbering[i][j] = renumber[i][oldrenum[j]];
        }
    }
}

void meshtracker::getrenumbering(std::shared_ptr<meshtracker> mt, std::vector<std::vector<int>>& renumbering)
{
    renumbering.resize(8);

    for (int i = 0; i < 8; i++)
    {
        renumbering[i].resize(elementrenumbering[i].size());
        if (mt != NULL && mt->elementrenumbering.size() != 0)
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

disjointregions* meshtracker::getdisjointregions(void)
{
    return &mydisjointregions;
}

void meshtracker::getindisjointregions(std::vector<std::vector<int>>& indisjregs)
{
    // Allocate the vectors:
    indisjregs.resize(8);
    for (int i = 0; i < 8; i++)
        indisjregs[i].resize(mydisjointregions.countelementsintype(i));
    
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

void meshtracker::print(void)
{
    for (int i = 0; i < elementrenumbering.size(); i++)
    {
        std::cout << std::endl << "Element type number " << i << ":" << std::endl << std::endl;
        for (int j = 0; j < elementrenumbering[i].size(); j++)
            std::cout << j << " --> " << elementrenumbering[i][j] << std::endl;
    }
    std::cout << std::endl;
}

