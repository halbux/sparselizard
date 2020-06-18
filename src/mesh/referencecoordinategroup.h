// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef REFERENCECOORDINATEGROUP_H
#define REFERENCECOORDINATEGROUP_H

#include <iostream>
#include <vector>
#include "myalgorithm.h"
#include "coordinategroup.h"

class referencecoordinategroup
{

    private:

        double noisethreshold = 1e-10;
        
        coordinategroup mycoordgroup;
        std::vector<int> myinputelems = {};
        std::vector<double> myinputrefcoords = {};
        
        // One entry per number of reference coordinates in an element:
        std::vector<std::vector<int>> myelems = {};
        std::vector<std::vector<int>> mycoordnums = {};
        std::vector<std::vector<double>> mykietaphis = {};
        
        
        // Current position status:
        int myrangebegin = 0, mynumrefcoords = 0;
        
        std::vector<double> mycurrefcoords = {};
        std::vector<int> mycurcoordnums = {};
        std::vector<int> mycurelems = {};
    
    public:
    
        referencecoordinategroup(void) {};
        referencecoordinategroup(std::vector<double>& coords);
        // Provide in 'elems' the element types and numbers {type0,elemnum0,type1,...}:
        referencecoordinategroup(std::vector<int>& elems, std::vector<double>& refcoords);
        
        // All disjoint regions should hold the same element type number:
        void evalat(std::vector<int> inputdisjregs);
        void evalat(int elemtypenum);
        
        void evalat(std::vector<int>& elems, std::vector<double>& kietaphis, std::vector<int>& coordnums);
        
        bool next(void);
        
        std::vector<double> getreferencecoordinates(void) { return mycurrefcoords; };
        std::vector<int> getcoordinatenumber(void) { return mycurcoordnums; };
        std::vector<int> getelements(void) { return mycurelems; };
        
};

#endif

