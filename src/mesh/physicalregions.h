// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef PHYSICALREGIONS_H
#define PHYSICALREGIONS_H

#include <iostream>
#include <vector>
#include "physicalregion.h"
#include "disjointregions.h"
#include <algorithm>
#include "myalgorithm.h"

class physicalregions
{

	private:

        std::vector<physicalregion> myphysicalregions;
        std::vector<int> myphysicalregionnumbers;
        
        disjointregions* mydisjointregions;
        
	public:
        
        physicalregions(disjointregions&);
                
        // Create a new physical region that is the union of all regions:
        int createunion(const std::vector<int> input);
        int createintersection(const std::vector<int> input);
		int createexclusion(int input, int toexclude);
        
        int getmaxphysicalregionnumber(void);
        
        // 'get' creates any non-existent physical region object:
        physicalregion* get(int physicalregionnumber);
        // 'getatindex' does not create non-existent physical region objects.
        physicalregion* getatindex(int physicalregionindex);
        // Get the number of physical regions:
        int count(void);
        // Get the total number of elements in all physical regions:
        int countelements(void);
        // Get all physical region numbers:
        std::vector<int> getallnumbers(void);
        // Get the physical region number of the physicalregionindex th physical region:
        int getnumber(int physicalregionindex);

};

#endif
