// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PETSCMESH_H
#define PETSCMESH_H

#include <string>
#include "element.h"
#include "nodes.h"
#include "elements.h"
#include "physicalregions.h"
#include "densematrix.h" 
#include "intdensematrix.h" 
#include "petscdmplex.h"
#include "petscviewer.h" 

class petscmesh
{
    private:

        DM mypetscmesh;
        
        int meshdim;
        
        int curvatureorder = 1;
        
        void reordernodes(int ourtypenum, std::vector<int>& toreorder);

    public:

        petscmesh(std::string filename);
        ~petscmesh(void);

        void extract(nodes& mynodes, elements& myelements, physicalregions& myphysicalregions, bool verbosity = false);
        
        // Write the petsc view to 'petscview.txt':
        void view(void);
};

#endif
