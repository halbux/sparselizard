// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPPARAMETER_H
#define OPPARAMETER_H

#include "operation.h"
#include "rawparameter.h"

class rawparameter;

class opparameter: public operation
{

    private:
        
        // Parameters are always reused
        
        int myrow;
        int mycolumn;
        
        std::shared_ptr<rawparameter> myparameter;
    
    public:
        
        opparameter(std::shared_ptr<rawparameter> input, int row, int col) { myparameter = input; myrow = row; mycolumn = col; };
        
        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        
        bool isparameter(void) { return true; };
        
        std::shared_ptr<rawparameter> getparameterpointer(void) { return myparameter; };
        
        int getselectedrow(void) { return myrow; };
        int getselectedcol(void) { return mycolumn; };
        
        bool isparameterincluded(std::vector<int> disjregs, rawparameter* rp);
        
        bool isharmonicone(std::vector<int> disjregs);
        
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        bool isvalueorientationdependent(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void print(void);

};

#endif
