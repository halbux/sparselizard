#include "disjointregionselector.h"


disjointregionselector::disjointregionselector(std::vector<int> disjointregionnumbers, std::vector<std::vector<int>> criteria)
{
    // Get the element type number in every disjoint region:
    std::vector<int> typenums(disjointregionnumbers.size());
    for (int i = 0; i < disjointregionnumbers.size(); i++)
        typenums[i] = (universe::mymesh->getdisjointregions())->getelementtypenumber(disjointregionnumbers[i]);
    criteria.push_back(typenums);
    
    std::vector<int> renumberingvector(disjointregionnumbers.size());
    std::iota(renumberingvector.begin(), renumberingvector.end(), 0);
    // Sort 'reorderingvector' according to all criteria.
    // The < operator is overloaded by a lambda function.
    std::sort(renumberingvector.begin(), renumberingvector.end(), [&](int elem1, int elem2)
    { 
        for (int i = 0; i < criteria.size(); i++)
        {
            if (criteria[i][elem1] < criteria[i][elem2])
                return true;
            if (criteria[i][elem1] > criteria[i][elem2])
                return false;
        }
        // For identical entries make a COHERENT decision for a stable sorting.
        return (elem1 < elem2);
    });
    
    // Sort all criteria accordingly:
    std::vector<std::vector<int>> criteriabackup = criteria;
    for (int i = 0; i < criteria.size(); i++)
    {
        for (int j = 0; j < criteria[i].size(); j++)
            criteria[i][j] = criteriabackup[i][renumberingvector[j]];
    }        
    
    // Split in groups:
    int currentgroup = -1;
    for (int i = 0; i < disjointregionnumbers.size(); i++)
    {
        bool newgroup = false;
        for (int j = 0; j < criteria.size(); j++)
            newgroup = newgroup || i == 0 || criteria[j][i] != criteria[j][i-1];
        
        if (newgroup)
        {
            groupeddisjointregions.push_back({disjointregionnumbers[renumberingvector[i]]});
            currentgroup++;
        }
        else
            groupeddisjointregions[currentgroup].push_back(disjointregionnumbers[renumberingvector[i]]);
    }
}
