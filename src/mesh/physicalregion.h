// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PHYSICALREGION_H
#define PHYSICALREGION_H

#include <iostream>
#include "disjointregions.h"
#include "physicalregions.h"
#include <vector>
#include <algorithm>
#include "element.h"

class physicalregions;

class physicalregion
{

    private:

        // The physical region can only hold elements of a single dimension (0D, 1D, 2D or 3D).
        int myelementdimension = -1;
        
        int myphysicalregionnumber;
        
        disjointregions* mydisjointregions;
        physicalregions* myphysicalregions;

        // 'includesdisjointregion[i]' is true if disjoint region i is in the physical region.
        std::vector<bool> includesdisjointregion = {};
        // List of all element numbers in the physical region.
        std::vector<std::vector<int>> elementlist = std::vector<std::vector<int>>(8, std::vector<int>(0));
        
    public:
        
        physicalregion(void) {};
        physicalregion(disjointregions&, physicalregions&, int physicalregionnumber, int elementdimension);
        
        int getnumber(void);
        // Add an element of uncurved type 'elementtypenumber' to the physical region:
        void addelement(int elementtypenumber, int elementnumber);
        
        int countelements(void);
        int getelementdimension(void);
        
        // Define the physical region in terms of the disjoint regions it contains:
        void definewithdisjointregions(int physregdim, std::vector<int> disjointregionlist, bool ismeshloading = false);

        // Get the definition of this physical region based on the disjoint regions it contains:
        std::vector<bool> getdefinition(void);

        // Get all disjoint regions of the max dimension:
        std::vector<int> getdisjointregions(void);
        // Get all disjoint regions of a given dimension (use -1 for all):
        std::vector<int> getdisjointregions(int dim);
        // Get all disjoint regions of a given element type:
        std::vector<int> getdisjointregionsoftype(int elementtypenumber);
        
        // 'renumberelements' updates the element numbers in 'elementlist' based on the input vector.
        // 'elementtypenumber' is the uncurved version of the actual element.
        void renumberelements(int elementtypenumber, std::vector<int>& elementrenumbering);
        
        // 'removeduplicatedelements' removes all duplicated elements in 'elementlist' (this call SORTS THE ELEMENTS).
        void removeduplicatedelements(void);
        
        // Get the elements in the physical region that have the region dimension:
        std::vector<std::vector<int>>* getelementlist(void);
        
        // Make a full copy of this object (linking objects used are the arguments):
        std::shared_ptr<physicalregion> copy(physicalregions* prs, disjointregions* drs);
        
};

#endif
