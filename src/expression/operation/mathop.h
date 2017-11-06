// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef MATHOP_H
#define MATHOP_H

#include <iostream>
#include "expression.h"
#include "integration.h"
#include "universe.h"
#include "vec.h"
#include "mat.h"

class expression;
class integration;
class mat;
class vec;

namespace mathop
{
    // Perform operations (union, intersection,...) on physical regions:
    int regionunion(const std::vector<int> physregs);
    int regionintersection(const std::vector<int> physregs);
    
    // Define the (unit norm) vector normal to a surface in 3D (or to a line in 2D). 
    expression normal(int surfphysreg);
    
    void setfundamentalfrequency(double f);
    void settime(double t);
    double gettime(void);
    
    // The time variable:
    expression t(void);
    
    expression dx(expression input);
    expression dy(expression input);
    expression dz(expression input);
    
    expression dt(expression input);
    expression dtdt(expression input);
    
    expression sin(expression input);
    expression cos(expression input);
    expression abs(expression input);
    expression sqrt(expression input);
    expression log10(expression input);
    expression pow(expression base, expression exponent);

    expression comp(int selectedcomp, expression input);
    expression compx(expression input);
    expression compy(expression input);
    expression compz(expression input);

    expression entry(int row, int col, expression input);

    expression transpose(expression input);
    expression inverse(expression input);
    expression determinant(expression input);

    expression grad(expression input);
    expression curl(expression input);
    
    integration integral(int physreg, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    integration integral(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    // For the multiharmonic resolution an extra integer is required 
    // to know with how many harmonics the coef multiplying the tf 
    // and/or dof should be approximated. Set 'numcoefharms' negative
    // and it will be as if you were calling the above non 
    // multiharmonic functions.
    integration integral(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    integration integral(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    
    expression dof(expression input, int physreg = -1);
    expression tf(expression input, int physreg = -1);
    
    
    
    // Define typically used arrays for convenience:
    expression array1x1(expression term11);
    expression array1x2(expression term11, expression term12);
    expression array1x3(expression term11, expression term12, expression term13);
    expression array2x1(expression term11, expression term21);
    expression array2x2(expression term11, expression term12, expression term21, expression term22);
    expression array2x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23);
    expression array3x1(expression term11, expression term21, expression term31);
    expression array3x2(expression term11, expression term12, expression term21, expression term22, expression term31, expression term32);
    expression array3x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23, expression term31, expression term32, expression term33);
    
    // Direct resolution:
    vec solve(mat A, vec b);/////////// CHECK AGAIN!!!!!!!!!!!!!!!!!!!
    
    
    
    ////////// PREDEFINED OPERATORS
    
    expression m2d(expression input);
    expression m3dn(expression input);
    expression m3ds(expression input);
    
    ////////// PREDEFINED FORMULATIONS
    
    expression predefinedelasticity(expression u, expression Eyoung, expression nupoisson);
    expression predefinedelectrostaticforce(expression gradtfu, expression gradv, expression epsilon);

};

#endif
