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
#include "gentools.h"

class physicalregion;

class physicalregions
{

    private:

        std::vector<std::shared_ptr<physicalregion>> myphysicalregions = {};
        std::vector<int> myphysicalregionnumbers = {};
        
        disjointregions* mydisjointregions = NULL;
        
    public:
        
        physicalregions(disjointregions&);
                
        // Create a new physical region that is the union of all regions:
        int createunion(std::vector<int> input, bool createifexisting = true);
        int createintersection(std::vector<int> input, int physregdim, bool createifexisting = true);
        int createunionofall(bool createifexisting = true);
        
        // Create a physical region from a list of disjoint regions:
        int createfromdisjointregionlist(int physregdim, std::vector<int> drs);
        
        // Define the physical regions based on the disjoint regions they contain:
        void definewithdisjointregions(void);
        
        int getmaxphysicalregionnumber(void);
        
        // 'get' creates any non-existent physical region object of requested dimension:
        physicalregion* get(int physicalregionnumber, int elementdimension = -1);
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
        
        // Find the physical region of given dimension and disjoint regions (-1 if not found):
        int find(int physregdim, std::vector<int> disjregsinphysreg);
        
        // Get the list of physical regions in which each element of a given type is.
        // 'addresses[i]' gives the first index in 'prs' where to find the physical
        // regions for the ith element. Length is numelems+1 and last entry gives prs.size().
        void inphysicalregions(int elementtypenumber, int totalnumelemsintype, std::vector<int>& addresses, std::vector<int>& prs);
        
        // Remove physical regions (do not call this yourself).
        void remove(std::vector<int> toremove, bool ispartofdisjregstructure);
        
        // Extract the positively renumbered physical regions:
        physicalregions extract(std::vector<int>& renumbering);
        
        // Give an error if any of the physical regions is not defined:
        void errorundefined(std::vector<int> physregs);
        // Give an error if not all physical regions have the same dimension:
        void errornotsamedim(std::vector<int> physregs);
        
        // Make a full copy of this object (linking objects used are the arguments):
        void copy(disjointregions* drs, physicalregions* target);
};

#endif
