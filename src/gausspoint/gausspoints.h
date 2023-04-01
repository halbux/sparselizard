// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef GAUSSPOINTS_H
#define GAUSSPOINTS_H

#include <iostream>
#include <vector>
#include <iomanip>

#include "logs.h"
#include "gppoint.h"
#include "gpline.h"
#include "gptriangle.h"
#include "gpquadrangle.h"
#include "gptetrahedron.h"
#include "gphexahedron.h"
#include "gpprism.h"
#include "gppyramid.h"


class gausspoints
{   
    private:
        
        int myelementtypenumber;
        
        std::vector<double> mycoordinates;
        std::vector<double> myweights;
        
    public:
        
        gausspoints(int elementtypenumber, int integrationorder);
        // Find based on the element type and the gauss points coordinates:
        gausspoints(int elementtypenumber, std::vector<double>& gpcoords);
        
        // 'getcoordinates' returns a vector containing the reference 
        // element coordinates of the Gauss points in the format 
        // [kigp1 etagp1 phigp1 kigp2 ...].
        std::vector<double> getcoordinates(void) { return mycoordinates; };
        std::vector<double> getweights(void) { return myweights; };
        
        int count(void) { return myweights.size(); };
        
        void print(void);
};

#endif
