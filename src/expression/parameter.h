// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include "operation.h"
#include "densematrix.h"
#include "disjointregions.h"
#include "universe.h"
#include "mesh.h"
#include "elementselector.h"
#include "parameterselectedregion.h"
#include "disjointregionselector.h"
#include "expression.h"
#include "field.h"
#include "vec.h"

class field;
class vec;
class operation;
class expression;
class parameterselectedregion;

class parameter
{

    private:

        int mynumrows;
        int mynumcols;

        // myexpressions[disjreg][i*mynumcols+j] stores the ith row, jth 
        // columns of the parameter expression defined on disjreg.
        std::vector<std::vector<std::shared_ptr<operation>>> myoperations;

        // Store a number associated to the operation on a given disjoint 
        // region. The operations on the disjoint regions on which the 
        // parameter has been defined with the same .set call have the same 
        // number because they are the same. This enables to interpolate on 
        // groups of disjoint regions that share the same operation number.
        int maxopnum = -1;
        std::vector<int> opnums;


        // Give an error if the parameter is undefined on at least one disj. reg.
        void errorifundefined(std::vector<int> disjregs);

        // Get the operation numbers for the requested disjoint regions.
        std::vector<int> getopnums(std::vector<int> disjregs);

    public:

        parameter(void);
        parameter(int numrows, int numcols);

        // A parameter cannot store an expression with a dof or a tf.
        // It can also only store expression arrays of a same dimension. 
        void set(int physreg, expression);

        std::shared_ptr<operation> get(int disjreg, int row, int col) { errorifundefined({disjreg}); return myoperations[disjreg][row*mynumcols+col]; };

        int countrows(void) { return mynumrows; };
        int countcolumns(void) { return mynumcols; };

        parameterselectedregion operator|(int physreg);

        std::vector<std::vector<densematrix>> interpolate(int row, int col, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int row, int col, int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        void simplify(int row, int col, int disjreg);

        void print(void);



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

        void write(int physreg, int numfftharms, std::string filename, int lagrangeorder = 1);
        void write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder = 1);

        void write(int physreg, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        void write(int physreg, expression meshdeform, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);


        // Defining the +, -, * and / operators:
        expression operator+(void);
        expression operator-(void);

        expression operator+(parameter&);
        expression operator-(parameter&);
        expression operator*(parameter&);
        expression operator/(parameter&);

        expression operator+(double);
        expression operator-(double);
        expression operator*(double);
        expression operator/(double);

};

// Define the left version of the operators based on the right one.
expression operator+(double, parameter&);
expression operator-(double, parameter&);
expression operator*(double, parameter&);
expression operator/(double, parameter&);

#endif
