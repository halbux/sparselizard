// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.



#ifndef SPANNINGTREE_H
#define SPANNINGTREE_H

#include <iostream>
#include <vector>
#include "mesh.h"
#include "elements.h"
#include <memory>

class spanningtree
{

    private:
        
        std::vector<int> startphysregs;
        
        elements* myelements;
        disjointregions* mydisjointregions;
        
        // Vector whose ith entry tells whether the edge number i is in the tree or not:
        std::shared_ptr<bool> isedgeintree = NULL;
        // The pointer to which the shared pointer points:
        bool* isedgeintreeptr;
        
        // Number of edges that are in the tree:
        int numberofedgesintree = 0;
        
        // True at index i if disjoint (edge) region i has priority in the tree construction:
        std::vector<bool> isprioritydisjointregion;
        
        
        // TEMPORARY CONTAINERS used during the tree construction.
        int numberofsubtrees = 0;
        // Entry i of this vector gives the subtree number in which edge i is:
        std::vector<int> insubtree;
        // List all edges in every subtree:
        std::vector<std::vector<int>> edgesinsubtree;
        // This vector gives true at entry i if subtree i has been added to the tree:
        std::vector<bool> issubtreeintree;
        // Entry i is true if node i is in tree:
        std::vector<bool> isnodeintree;
        
        
        
        // First phase of the tree creation. All (unconnected) subtrees 
        // are created on the priority disjoint edge regions.
        // Subtrees are not allowed to share nodes or edges with each other.
        void growsubtrees(void);
        // Grow the subtree that has edges only on the priority edge disjoint regions and 
        // starting at node 'nodenumber'. Give it subtree number 'subtreenumber'. 
        // This can only be called if at least one edge can be added to the subtree.
        void growsubtree(int nodenumber, int subtreenumber);
        
        // Create the final tree by connecting all subtrees together:
        void connectsubtrees(void);
        // Grow the tree starting at a given node.
        // The subtrees must have been defined before the call.
        void growtree(int nodenumber);
        
        void grow(void);
        
        
        int mymeshnumber = 0;
        // Synchronize with the hp-adapted mesh:
        void synchronize(void);
        // To avoid infinite recursive calls:
        bool issynchronizing = false;

    public:

        // Create a spanning tree on the whole mesh.
        // The tree is first created on the physical regions provided.
        spanningtree(std::vector<int> physregs);
        
        // Is the index th element of the disjoint region in the tree? Always false if not an edge.
        bool isintree(int index, int disjreg);
        
        // Count edges in a disjoint edge region of the tree (return 0 if not a disjoint edge region):
        int countedgesintree(int disjreg);
        
        // Count edges in the tree:
        int countedgesintree(void);
        // Get all edges in the tree (ordering does not follow tree):
        std::vector<int> getedgesintree(void);
        
        // Write to file to visualize the tree:
        void write(std::string filename);
    
};

#endif
 
