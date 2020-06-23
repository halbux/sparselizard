// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object keeps track of the element numbering when the mesh is p-adapted.

#ifndef PTRACKER_H
#define PTRACKER_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "disjointregions.h"

class ptracker
{

    private:
        
        disjointregions mydisjointregions;
        
        // The elements of type i can be renumbered as new number = elementrenumbering[i][old number].
        std::vector<std::vector<int>> elementrenumbering = {};
    
    public:

        ptracker(std::vector<int> numelemspertype);
        
        // This updates the disjoint region object:
        void updatedisjointregions(disjointregions* input);
        
        // Update the stored renumbering according to the subsequent renumbering step provided as argument.
        void updaterenumbering(std::vector<std::vector<int>>& renumber);
        
        // The element number in the mesh tracked by this mesh tracker can be renumbered to the one
        // in the mesh tracked by the 'mt' tracker by using mtelemnum = renumbering[thiselemnum].
        // In case 'mt' is NULL it is considered to have identity renumbering (no number change).
        void getrenumbering(std::shared_ptr<ptracker> mt, std::vector<std::vector<int>>& renumbering);
        
        disjointregions* getdisjointregions(void);
        
        // Provides a vector indisjregs[i][e] giving the disjoint region in which element e of type i is (in the mesh state of this mesh tracker).
        void getindisjointregions(std::vector<std::vector<int>>& indisjregs);
        
        void print(void);
    
};

#endif

