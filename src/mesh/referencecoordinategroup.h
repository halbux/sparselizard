// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef REFERENCECOORDINATEGROUP_H
#define REFERENCECOORDINATEGROUP_H

#include <iostream>
#include <vector>
#include "coordinategroup.h"

class referencecoordinategroup
{

    private:

        double noisethreshold = 1e-8;
        
        int mynumcoords = 0;
        coordinategroup mycoordgroup;
        
        std::vector<int> mydisjregs = {};
        
        std::vector<int> myelems = {};
        std::vector<double> mykietaphis = {};
        
        // Working range begin and end:
        int myrangebegin = 0;
        int myrangelength = 0;
        
        // To find back the original coordinate numbers:
        std::vector<int> myunorderingvector = {};
        
        std::vector<double> mycurrefcoords = {};
        std::vector<int> mycurcoordnums = {};
        std::vector<int> mycurelems = {};
    
    public:
    
        referencecoordinategroup(void) {};
        referencecoordinategroup(std::vector<double>& coords);
        
        void evalat(std::vector<int> inputdisjregs);
        
        void group(void);
        
        bool next(void);
        
        std::vector<double> getreferencecoordinates(void) { return mycurrefcoords; };
        std::vector<int> getcoordinatenumber(void) { return mycurcoordnums; };
        std::vector<int> getelements(void) { return mycurelems; };
        
};

#endif

