// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef REGIONDEFINER_H
#define REGIONDEFINER_H

#include <iostream>
#include <vector>
#include "element.h"
#include "elements.h"
#include "nodes.h"
#include "physicalregions.h"

class regiondefiner
{
    private:
    
        double roundoffnoise = 1e-10;

        nodes* mynodes;
        elements* myelements;
        physicalregions* myphysicalregions;

        // mypriority[i] gives {operationtype,operationindex} for operation of priority i.
        std::vector<std::vector<int>> mypriority = {};


        // SKIN (operation type 0)
        //
        // List of new skin physical regions:
        std::vector<int> skins = {};
        // List of physical regions for which the skin is requested:
        std::vector<int> toskin = {};

        // BOX SELECTION (operation type 1)
        //
        // List of new box-selected physical regions:
        std::vector<int> boxed = {};
        // List of physical regions for which the box selection is requested:
        std::vector<int> tobox = {};
        // List of element dimensions to select:
        std::vector<int> boxelemdims = {};
        // List of the x, y and z box limits {x1,x2,y1,y2,z1,z2}:
        std::vector<std::vector<double>> boxlimits = {};

        // SPHERE SELECTION (operation type 2)
        //
        // List of new sphere-selected physical regions:
        std::vector<int> sphered = {};
        // List of physical regions for which the sphere selection is requested:
        std::vector<int> tosphere = {};
        // List of element dimensions to select:
        std::vector<int> sphereelemdims = {};
        // List of the x, y and z sphere centers {xc,yc,zc}:
        std::vector<std::vector<double>> spherecenters = {};
        // List of the sphere radii:
        std::vector<double> sphereradii = {};


        void defineskinregion(int regnum);
        void defineboxregion(int regnum);
        void definesphereregion(int regnum);

    public:

        regiondefiner(nodes& inputnodes, elements& inputelems, physicalregions& inputphysregs);

        void regionskin(int newphysreg, int physregtoskin);
        void boxselection(int newphysreg, int selecteddim, std::vector<double> boxlimit, int physregtobox);
        void sphereselection(int newphysreg, int selecteddim, std::vector<double> centercoords, double radius, int physregtosphere);


        void defineregions(void);
};

#endif
