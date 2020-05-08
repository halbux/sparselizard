// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef LAGRANGEFORMFUNCTION_H
#define LAGRANGEFORMFUNCTION_H

#include <iostream>
#include <vector>
#include "densematrix.h"
#include "polynomial.h"
#include "element.h"
#include <iomanip>
#include "universe.h"

#include "lagrangepoint.h"
#include "lagrangeline.h"
#include "lagrangetriangle.h"
#include "lagrangequadrangle.h"
#include "lagrangetetrahedron.h"
#include "lagrangehexahedron.h"
#include "lagrangeprism.h"
#include "lagrangepyramid.h"


class lagrangeformfunction
{
    private:
        
        int myorder;
        int myelementtypenumber;
        std::vector<double> myevaluationpoints;
        
        std::vector<double> mynodecoordinates = {};
        std::vector<polynomial> myformfunctionpolynomials = {};
        
        // In the following dense matrices there is a row per Lagrange 
        // form function and a column per evaluation point.
        // evaluated[i] stores the value of the Lagrange form functions:
        // - without  derivative for i = 0
        // - with ki  derivative for i = 1
        // - with eta derivative for i = 2
        // - with phi derivative for i = 3
        std::vector<densematrix> evaluated = std::vector<densematrix>(4);
        
        // Make the Lagrange polynomials/coordinates available: 
        void preparepoly(void);
        void preparecoords(void);
        
    public:
        
        lagrangeformfunction(void) {};
        
        // The constructor takes the element type number, the form function
        // order and the coordinates of the evaluation points as input. 
        // They are set once and for all for a given object. 
        // Format is [ki1 eta1 phi1 ki2 eta2 phi2 ki3 ...].
        lagrangeformfunction(int elementtypenumber, int order, const std::vector<double> evaluationpoints);
        
        // getderivative gives the value of the Lagrange form functions:
        // - without  derivative for whichderivative = 0
        // - with ki  derivative for whichderivative = 1
        // - with eta derivative for whichderivative = 2
        // - with phi derivative for whichderivative = 3
        //
        // It outputs a copy of evaluated[whichderivative] if available.
        densematrix getderivative(int whichderivative);
        
        std::vector<double> getnodecoordinates(void);
        std::vector<polynomial> getformfunctionpolynomials(void);
        
        void print(void);
};

#endif
