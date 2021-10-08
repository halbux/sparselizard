// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include "universe.h"
#include "rawparameter.h"
#include "parameterselectedregion.h"
#include "expression.h"
#include "field.h"
#include "vec.h"

class field;
class vec;
class expression;
class parameterselectedregion;

class parameter
{

    private:

        // The actual parameter:
        std::shared_ptr<rawparameter> rawparamptr = NULL;
    
    public:

        parameter(void);
        parameter(int numrows, int numcols);

        int countrows(void);
        int countcolumns(void);

        parameterselectedregion operator|(int physreg);

        void print(void);

        std::shared_ptr<rawparameter> getpointer(void) { return rawparamptr; };


        vec atbarycenter(int physreg, field onefield);

        std::vector<double> max(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        
        void interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        void interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        
        std::vector<double> interpolate(int physreg, const std::vector<double> xyzcoord);
        std::vector<double> interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord);

        double integrate(int physreg, expression meshdeform, int integrationorder);
        double integrate(int physreg, int integrationorder);

        void write(int physreg, int numfftharms, std::string filename, int lagrangeorder);
        void write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder);

        void write(int physreg, std::string filename, int lagrangeorder, int numtimesteps = -1);
        void write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps = -1);


        // Defining the +, -, * and / operators:
        expression operator+(void);
        expression operator-(void);

        expression operator+(parameter);
        expression operator-(parameter);
        expression operator*(parameter);
        expression operator/(parameter);

        expression operator+(double);
        expression operator-(double);
        expression operator*(double);
        expression operator/(double);

};

// Define the left version of the operators based on the right one.
expression operator+(double, parameter);
expression operator-(double, parameter);
expression operator*(double, parameter);
expression operator/(double, parameter);

#endif
