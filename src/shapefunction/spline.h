// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// This objects holds a cubic natural spline interpolation of a dataset.


#ifndef SPLINE_H
#define SPLINE_H

#include <iostream>
#include <vector>
#include <string>
#include "myalgorithm.h"
#include "densematrix.h"
#include "mathop.h"
#include "mat.h"
#include "vec.h"

class spline
{

    private:

        double noisethreshold = 1e-10;

        double xmin, xmax;
        // The x and y axis data (sorted ascendingly):
        densematrix myx, myy;
        // The spline parameters:
        densematrix mya, myb;
    
    public:
    
        spline(void) {};
        spline(std::string filename, char delimiter = '\n');
        spline(std::vector<double>& xin, std::vector<double>& yin);
        
        void set(std::vector<double>& xin, std::vector<double>& yin);
        
        densematrix evalat(densematrix input);
        
        void write(std::string filename, int numsplits, char delimiter = '\n');
    
};

#endif

