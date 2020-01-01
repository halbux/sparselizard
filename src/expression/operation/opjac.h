// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPJAC_H
#define OPJAC_H

#include "operation.h"

class opjac: public operation
{

    private:
        
        int myrow, mycol;
    
    public:
        
        // Define as jac(m,n):
        opjac(int m, int n) { myrow = m; mycol = n; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);
        
        void print(void);
        
};

#endif
