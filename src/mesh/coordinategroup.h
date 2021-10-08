// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef COORDINATEGROUP_H
#define COORDINATEGROUP_H

#include <iostream>
#include <vector>
#include <memory>

class coordinategroup
{

    private:

        int mynumcoords = 0;
        
        // Bounds gives {xmin, xmax, ymin, ymax, zmin, zmax}:
        std::vector<double> bounds = {};
        // Delta gives {deltax, deltay, deltaz} between two slice tics:
        std::vector<double> delta = {};
        // Number of x, y and z slices:
        std::vector<int> numslices = {};
        
        std::vector<double> noisethreshold = {};
    
        // First index in 'mygroups' of each group (last element is number of coordinates):
        std::vector<int> mygroupads = {};
        // Coordinate indexes in each group:
        std::shared_ptr<int> mygroups;
        // Coordinates in each group:
        std::shared_ptr<double> mygroupcoords;
        
        int selx1, selx2, sely1, sely2, selz1, selz2;
        int curselx, cursely, curselz, curg;
    
    public:
    
        coordinategroup(void) {};
        coordinategroup(std::vector<double>& coords);

        int countcoordinates(void);
        
        void select(double x, double y, double z, double radius);
        // Move to next group. Return false if none.
        bool next(void);
        
        int countgroupcoordinates(void);
        int* getgroupindexes(void);
        double* getgroupcoordinates(void);
            
};

#endif

