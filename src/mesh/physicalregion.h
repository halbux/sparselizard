// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef PHYSICALREGION_H
#define PHYSICALREGION_H

#include <iostream>
#include "disjointregions.h"
#include <vector>
#include <algorithm>
#include "element.h"

class physicalregion
{

	private:

        // The physical region can only hold elements of a single dimension (0D, 1D, 2D or 3D).
        int myelementdimension = -1;
        
        int myphysicalregionnumber;
        int myphysicalregionindex;
        
        disjointregions* mydisjointregions;

        // 'includesdisjointregion[i]' is true if disjoint region i is in the physical region.
        std::vector<bool> includesdisjointregion;
        // List of all element numbers for every type in every physical region defined in the loaded mesh file.
        std::vector<std::vector<int>> elementlist = std::vector<std::vector<int>>(8, std::vector<int>(0));
        
	public:
        
        physicalregion(disjointregions&, int physicalregionnumber, int physicalregionindex);
        
        int getnumber(void);
        // Add an element of uncurved type 'elementtypenumber' to the physical region:
        void addelement(int elementtypenumber, int elementnumber);
        
        int countelements(void);
        int getelementdimension(void);
        
        // 'definewithdisjointregions' defines the physical region in terms
        // of the disjoint regions it contains. Note: calling it clears 'elementlist'.
        void definewithdisjointregions(void);  
        // Set manually the list of disjoint regions in the physical region.
        // Duplicates are removed. 'myelementdimension' is updated to the max dim of the disj regs.
        void setdisjointregions(std::vector<int> disjointregionlist);

        // Get all disjoint regions of the max dimension:
        std::vector<int> getdisjointregions(void);
        // Get all disjoint regions of a given dimension (use -1 for all):
        std::vector<int> getdisjointregions(int dim);
        
        // 'renumberelements' updates the element numbers in 'elementlist' based on the input vector.
        // 'elementtypenumber' is the uncurved version of the actual element.
        void renumberelements(int elementtypenumber, std::vector<int>& elementrenumbering);
        
        // 'removeduplicatedelements' removes all duplicated elements in 'elementlist'. 
        // As a side effect all elements are sorted in a every physical region.
        void removeduplicatedelements(void);
        
        std::vector<std::vector<int>>* getelementlist(void);
        
};

#endif
