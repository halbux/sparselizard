// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef SPANNINGTREE_H
#define SPANNINGTREE_H

#include <iostream>
#include <vector>
#include <memory>
#include "rawspanningtree.h"

class rawspanningtree;

class spanningtree
{
    private:
        
        std::shared_ptr<rawspanningtree> myrawspantree = NULL;

    public:

        // Create a spanning tree on the whole mesh. The tree flows from
        // the highest (first) to the lowest (last) priority region.
        spanningtree(std::vector<int> physregs);
        
        // Count the number of edges in the tree:
        int countedgesintree(void);
        // Get all edges in the tree (ordering does not follow tree):
        std::vector<int> getedgesintree(void);
        
        // Get the raw spanning tree pointer:
        std::shared_ptr<rawspanningtree> getpointer(void);
        
        // Write to file to visualize the tree:
        void write(std::string filename);
};

#endif
 
