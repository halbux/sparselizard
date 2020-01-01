// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// This object holds a cubic natural spline interpolation of a dataset.


#ifndef SPLINE_H
#define SPLINE_H

#include <iostream>
#include <vector>
#include <string>
#include "densematrix.h"

class spline
{

    private:

        double noisethreshold = 1e-10;

        double xmin, xmax;
        // The x and y axis data (sorted ascendingly):
        densematrix myx, myy;
        // The spline parameters:
        densematrix mya, myb;
        
        int derivativeorder = 0;
    
    public:
    
        spline(void) {};
        spline(std::string filename, char delimiter = '\n');
        spline(std::vector<double> xin, std::vector<double> yin);
        
        void set(std::vector<double>& xin, std::vector<double>& yin);
        
        spline getderivative(void);
        
        double getxmin(void) { return xmin; };
        double getxmax(void) { return xmax; };
        
        double evalat(double input);
        std::vector<double> evalat(std::vector<double> input);
        densematrix evalat(densematrix input);
        
        void write(std::string filename, int numsplits, char delimiter = '\n');
    
};

#endif

