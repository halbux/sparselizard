// sparselizard - Copyright (C) see copyright file.
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
        
        // EXCLUSION (operation type 3)
        //
        // List of new exclusion physical regions:
        std::vector<int> excluded = {};
        // List of physical regions from which an exclusion is requested:
        std::vector<int> toexcludefrom = {};
        // List of physical regions to exclude:
        std::vector<std::vector<int>> toexclude = {};

        // LAYER SELECTION (operation type 4)
        //
        // List of new layer-selected physical regions:
        std::vector<int> layered = {};
        // List of physical regions from which the layer selection is requested:
        std::vector<int> tolayer = {};
        // List of physical regions where the layer selection growth starts:
        std::vector<int> growthstart = {};
        // List of number of layers to select:
        std::vector<int> numlayers = {}; 
        
        // ANY NODE SELECTION (operation type 5)
        //
        // List of new node regions:
        std::vector<int> anynoded = {};
        // List of physical regions from which any node is requested:
        std::vector<int> toanynode = {};
        

        void defineskinregion(int regnum);
        void defineboxregion(int regnum);
        void definesphereregion(int regnum);
        void defineexclusionregion(int regnum);
        void definelayerregion(int regnum);
        void defineanynoderegion(int regnum);

    public:

        regiondefiner(void) {};
        regiondefiner(nodes& inputnodes, elements& inputelems, physicalregions& inputphysregs);

        void regionskin(int newphysreg, int physregtoskin);
        void regionbox(int newphysreg, int selecteddim, std::vector<double> boxlimit, int physregtobox);
        void regionsphere(int newphysreg, int selecteddim, std::vector<double> centercoords, double radius, int physregtosphere);
        void regionexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude);
        void regionlayer(int newphysreg, int physregtoselectfrom, int physregtostartgrowth, int nl);
        void regionanynode(int newphysreg, int physregtoselectfrom);

        bool isanyregiondefined(void);
        bool isanycoordinatedependentregiondefined(void);

        void defineregions(void);
        
        // Make a full copy of this object (linking objects used are the arguments):
        regiondefiner copy(nodes* nds, elements* els, physicalregions* prs);
};

#endif
