// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object holds a geometrical shape.

#ifndef RAWEXTRUSION_H
#define RAWEXTRUSION_H

#include <iostream>
#include <vector>

#include "expression.h"

#include "rawshape.h"
#include "rawquadrangle.h"
#include "rawtriangle.h"
#include "rawline.h"
#include "rawpoint.h"
#include "geotools.h"

class rawextrusion: public rawshape
{   

    private:

        int myphysicalregion = -1;
        
        // Son shapes:
        std::vector<std::shared_ptr<rawshape>> sons = {};

        // Coordinates of the nodes in the mesh:
        std::vector<double> mycoords = {};
        // Elements in the mesh:
        std::vector<std::vector<int>> myelems = std::vector<std::vector<int>>(8, std::vector<int>(0));


        // Number of node layers in the extrusion:
        int mynumlayers;

        // Extrusion length:
        double myheight;
        
        // Extrusion direction:
        std::vector<double> myextrudedirection;

        // Unextruded rawshape:
        std::shared_ptr<rawshape> mybaseshape;
        
    public:

        rawextrusion(void) {};

        rawextrusion(int physreg, std::shared_ptr<rawshape> innerrawshape, double height, int numlayers, std::vector<double> extrudedirection);

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

        void mesh(void);

};

#endif
