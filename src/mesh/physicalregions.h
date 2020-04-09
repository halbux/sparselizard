// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PHYSICALREGIONS_H
#define PHYSICALREGIONS_H

#include <iostream>
#include <vector>
#include <memory>
#include "physicalregion.h"
#include "disjointregions.h"
#include <algorithm>
#include "myalgorithm.h"

class physicalregion;

class physicalregions
{

    private:

        std::vector<std::shared_ptr<physicalregion>> myphysicalregions;
        std::vector<int> myphysicalregionnumbers;
        
        disjointregions* mydisjointregions;
        
    public:
        
        physicalregions(disjointregions&);
                
        // Create a new physical region that is the union of all regions:
        int createunion(const std::vector<int> input);
        int createintersection(const std::vector<int> input);
        int createexclusion(int input, int toexclude);
        int createunionofall(void);
        
        int getmaxphysicalregionnumber(void);
        
        // 'get' creates any non-existent physical region object:
        physicalregion* get(int physicalregionnumber);
        // 'getatindex' does not create non-existent physical region objects.
        physicalregion* getatindex(int physicalregionindex);
        // Get the number of physical regions of a given dimension (use -1 for all):
        int count(int dim = -1);
        // Get the total number of elements in all physical regions:
        int countelements(void);
        // Get all physical region numbers of a given dimension (use -1 for all):
        std::vector<int> getallnumbers(int dim = -1);
        // Get the physical region number of the physicalregionindex th physical region:
        int getnumber(int physicalregionindex);
        // Get the index of the physical region number (-1 if undefined):
        int getindex(int physicalregionnumber);
        
        // Remove physical regions (physical regions required for the compressed mesh structure should NOT be removed).
        void remove(std::vector<int> toremove, bool ispartofdisjregstructure);
        
        // Give an error if any of the physical regions is not defined:
        void errorundefined(std::vector<int> physregs);
        
        // Clear the object (pointers to other objects are not cleared).
        void clear(void);
};

#endif
