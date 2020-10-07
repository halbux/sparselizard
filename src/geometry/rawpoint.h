// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object holds a geometrical shape.

#ifndef RAWPOINT_H
#define RAWPOINT_H

#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>

#include "rawshape.h"
#include "rawextrusion.h"

class rawpoint: public rawshape
{   

    private:
        
        int myphysicalregion = -1;

        // Son shapes:
        std::vector<std::shared_ptr<rawshape>> sons = {};

        // Coordinates of the nodes in the mesh:
        std::vector<double> mycoords = {};
        // Elements in the mesh:
        std::vector<std::vector<int>> myelems = std::vector<std::vector<int>>(8, std::vector<int>(0));

    public:
 
        rawpoint(void) {};

        rawpoint(int physreg, std::vector<double> coords);
      
        // Provide to this constructor the coordinates of all mesh nodes and elements in the point cloud:
        rawpoint(int physreg, std::vector<double>& allcoords, std::vector<std::vector<int>>& allelems);
        
        std::shared_ptr<rawshape> extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection);

        std::shared_ptr<rawshape> duplicate(void);

        void setphysicalregion(int physreg);
    
        int getdimension(void);

        std::string getname(void);

        std::vector<std::shared_ptr<rawshape>> getsons(void);

        // Get subshapes (sons are included):
        std::vector<std::shared_ptr<rawshape>> getsubshapes(void);
        void setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes);

        // Get the mesh info of the shape:
        int getphysicalregion(void);
        std::vector<double>* getcoords(void); 
        std::vector<std::vector<int>>* getelems(void);


        std::shared_ptr<rawshape> getpointer(void);
};

#endif
