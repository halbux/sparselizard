// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef COORDINATEGROUP_H
#define COORDINATEGROUP_H

#include <iostream>
#include <vector>
#include "universe.h"

class coordinategroup
{

    private:

        int mynumcoords = 0;

        // Target N elements per block on a structured grid:
        int N = 100;
        
        // Bounds gives {xmin, xmax, ymin, ymax, zmin, zmax}:
        std::vector<double> bounds = {};
        // Delta gives {deltax, deltay, deltaz} between two slice tics:
        std::vector<double> delta = {};
        // Number of x, y and z slices:
        std::vector<int> numslices = {};
        
        std::vector<double> noisethreshold = {};
    
        // Element indexes in each group:
        std::vector<std::vector<std::vector<std::vector<int>>>> mygroups = {};
        std::vector<std::vector<std::vector<std::vector<double>>>> mygroupcoords = {};
        

        // Selected coordinate to locate in the boxes:        
        double xselect, yselect, zselect, myradius;
        
        std::vector<int> selectedgroups = {};
    
    public:

        coordinategroup(std::vector<double>& coords);

        void select(double x, double y, double z, double maxelemsize);

        int countgroups(void);
                
        std::vector<int>* getgroupindexes(int groupnum);
        std::vector<double>* getgroupcoordinates(int groupnum);
            
};

#endif

