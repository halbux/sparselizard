// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef RAWPARAMETER_H
#define RAWPARAMETER_H

#include <iostream>
#include "operation.h"
#include "densemat.h"
#include "universe.h"
#include "elementselector.h"
#include "expression.h"

class operation;

class rawparameter
{

    private:

        int mynumrows;
        int mynumcols;

        // myexpressions[disjreg][i*mynumcols+j] stores the ith row, jth 
        // columns of the parameter expression defined on disjreg.
        std::vector<std::vector<std::shared_ptr<operation>>> myoperations = {};

        // Store a number associated to the operation on a given disjoint 
        // region. The operations on the disjoint regions on which the 
        // parameter has been defined with the same .set call have the same 
        // number because they are the same. This enables to interpolate on 
        // groups of disjoint regions that share the same operation number.
        int maxopnum = -1;
        std::vector<int> opnums = {};
        
        
        int mymeshnumber = 0;
        
        // Track the calls to 'set'.
        std::vector<std::pair<int, expression>> mystructuretracker = {};
        
        // Synchronize with the hp-adapted mesh:
        void synchronize(void);
        // To avoid infinite recursive calls:
        bool issynchronizing = false;


        // Give an error if the parameter is undefined on at least one disj. reg.
        void errorifundefined(std::vector<int> disjregs);

        // Get the operation numbers for the requested disjoint regions.
        std::vector<int> getopnums(std::vector<int> disjregs);

    public:

        rawparameter(int numrows, int numcols);

        // A parameter cannot store an expression with a dof or a tf.
        // It can also only store expression arrays of a same dimension. 
        void set(int physreg, expression input);
        
        bool isdefined(int disjreg);
        
        std::shared_ptr<operation> get(int disjreg, int row, int col);

        int countrows(void);
        int countcolumns(void);

        std::vector<std::vector<densemat>> interpolate(int row, int col, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int row, int col, int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        void simplify(int row, int col, int disjreg);

        void print(void);

};

#endif
