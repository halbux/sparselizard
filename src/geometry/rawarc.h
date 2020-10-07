// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object holds a geometrical shape.

#ifndef RAWARC_H
#define RAWARC_H

#include <iostream>
#include <vector>
#include <memory>

#include "geotools.h"
#include "rawshape.h"
#include "rawpoint.h"
#include "rawextrusion.h"

class rawarc: public rawshape
{   

    private:
        
        int myphysicalregion = -1;
        
        int mynummeshpoints;

        // Son shapes:
        std::vector<std::shared_ptr<rawshape>> sons = {};
        // Arc center point:
        std::shared_ptr<rawshape> mycenterpoint;

        // Coordinates of the nodes in the mesh:
        std::vector<double> mycoords = {};
        // Elements in the mesh:
        std::vector<std::vector<int>> myelems = std::vector<std::vector<int>>(8, std::vector<int>(0));

    public:

        rawarc(void) {};
    
        // Provide to this constructor the point rawshapes at the two ends of the arc and at its center.
        // Order is start point, end point, center point.
        rawarc(int physreg, std::vector<std::shared_ptr<rawshape>> inputpoints, int nummeshpoints);

        std::shared_ptr<rawshape> extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection);

        std::shared_ptr<rawshape> duplicate(void);

        // Flip the direction:
        void flip(void);

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

        void mesh(void);

};

#endif
