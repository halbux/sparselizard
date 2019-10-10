// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// This objects holds a cubic natural spline interpolation of a dataset.


#ifndef SPLINES_H
#define SPLINES_H

#include <iostream>
#include <vector>
#include <string>
#include "myalgorithm.h"
#include "densematrix.h"
#include "mathop.h"
#include "mat.h"
#include "vec.h"

class splines
{

    private:

        // The x and y axis data (sorted ascendingly):
        densematrix myx, myy;
        // The spline parameters:
        densematrix myy0,myy1,mya,myb;
    
    public:

        splines(std::vector<double>& xin, std::vector<double>& yin);
        
        densematrix evalat(densematrix input);
        
        void write(std::string filename, int numevalsperspline);
    
};

#endif

