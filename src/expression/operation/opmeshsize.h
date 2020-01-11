// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPMESHSIZE_H
#define OPMESHSIZE_H

#include "operation.h"

class opmeshsize: public operation
{

    private:
        
        bool reuse = false;
        
        int myintegrationorder;
        
    public:
        
        opmeshsize(int integrationorder) { myintegrationorder = integrationorder; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
