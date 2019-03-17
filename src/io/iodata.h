// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef IODATA_H
#define IODATA_H

#include <iostream>
#include <vector>
#include "densematrix.h"

class iodata
{   
    
    private:
    
        // Interpolation order for the field data and for the geometry:
        int myinterpolorder;
        int mygeointerpolorder;
        
        // Scalar or vector data (3 components):
        bool isscalardata;
        
        // IN CASE OF MULTIPLE TIMESTEPS THE DATA BLOCKS BELOW ARE 
        // THE COLUMN-WISE CONCATENATION OF SINGLE-TIMESTEP DATA.
        //
        // Time vals in case of data at multiple timesteps.
        std::vector<double> mytimevals = {};
        
        // Coordinates of the nodes at which to write the data.
        // mycoords[xyz][i] gives the x, y or z coordinates for all elements of type i.
        // Every row corresponds to a given element. Every column corresponds to a node in the element.
        std::vector<std::vector<std::vector<densematrix>>> mycoords;
        
        // mydata[comp][i] gives the data at component 'comp' for all elements of type i.
        // The data for scalars is at component 0.
        std::vector<std::vector<std::vector<densematrix>>> mydata;		
        
        // Combine the data in every element type:
        void combine(void);
    
    public:
    
        iodata(int interpolorder, int geointerpolorder, bool isitscalardata, std::vector<double> timevals);
        
        bool isscalar(void);
        int getinterpolorder(void);
        int getgeointerpolorder(void);
        
        // Get time tag:
        std::vector<double> gettimetags(void);
        
        // ALWAYS CHECK THAT THERE IS DATA FOR A GIVEN ELEMENT TYPE BEFORE REQUESTING IT:
        bool ispopulated(int elemtypenum);
        
        // Rows correspond to the elements. Columns correspond to the element nodes.
        void addcoordinates(int elemtypenum, densematrix xcoords, densematrix ycoords, densematrix zcoords);
        // Add all components vals[comp] of the data:
        void adddata(int elemtypenum, std::vector<densematrix> vals);
        
        // Get the x, y and z coordinates of the nodes in all elements of a given type:
        std::vector<densematrix> getcoordinates(int elemtypenum);
        // Get all components of the data at the nodes in all elements of a given type:
        std::vector<densematrix> getdata(int elemtypenum);
        
};

#endif
