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
        int createunion(std::vector<int> input, bool createifexisting = true);
        int createintersection(std::vector<int> input, bool createifexisting = true);
        int createunionofall(bool createifexisting = true);
        
        // Create a physical region from a list of disjoint regions of same dimension:
        int createfromdisjointregionlist(std::vector<int> drs);
        
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
        
        // Get the number of the physical region made exactly of the argument disjoint regions (-1 if not found):
        int find(std::vector<int>& disjregsinphysreg);
        
        // Get the list of physical regions in which each element of a given type is.
        // 'addresses[i]' gives the first index in 'prs' where to find the physical
        // regions for the ith element. Length is numelems+1 and last entry gives prs.size().
        void inphysicalregions(int elementtypenumber, int totalnumelemsintype, std::vector<int>& addresses, std::vector<int>& prs);
        
        // Remove physical regions (do not call this yourself).
        void remove(std::vector<int> toremove, bool ispartofdisjregstructure);
        
        // Give an error if any of the physical regions is not defined:
        void errorundefined(std::vector<int> physregs);
        // Give an error if not all physical regions have the same dimension:
        void errornotsamedim(std::vector<int> physregs);
        
        // Make a full copy of this object (linking objects used are the arguments):
        void copy(disjointregions* drs, physicalregions* target);
};

#endif
