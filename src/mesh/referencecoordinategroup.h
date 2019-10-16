// sparselizard - Copyright (C) 2017- A. Halbach
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

        double noisethreshold = 1e-8;
        
        coordinategroup mycoordgroup;
        
        std::vector<int> mydisjregs = {};
        
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
        
        // All disjoint regions should hold the same element type number:
        void evalat(std::vector<int> inputdisjregs);
        
        bool next(void);
        
        std::vector<double> getreferencecoordinates(void) { return mycurrefcoords; };
        std::vector<int> getcoordinatenumber(void) { return mycurcoordnums; };
        std::vector<int> getelements(void) { return mycurelems; };
        
};

#endif

