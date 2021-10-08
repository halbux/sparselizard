// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef DISJOINTREGIONSELECTOR_H
#define DISJOINTREGIONSELECTOR_H

#include <iostream>
#include "disjointregions.h"
#include "universe.h"
#include <vector>
#include <numeric>

class disjointregionselector
{
    private:
    
        // 'groupeddisjointregions[i]' holds the ith group.
        std::vector<std::vector<int>> groupeddisjointregions = {};
        
    public:
    
        // Group by same element type number and same value for all criteria[i]:
        disjointregionselector(std::vector<int> disjointregionnumbers, std::vector<std::vector<int>> criteria);
        
        int countgroups(void) { return groupeddisjointregions.size(); };
        std::vector<int> getgroup(int groupnumber) { return groupeddisjointregions[groupnumber]; };
        
};

#endif
