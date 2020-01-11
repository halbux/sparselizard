// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPFIELDORDER_H
#define OPFIELDORDER_H

#include "operation.h"

class opfieldorder: public operation
{

    private:
        
        bool reuse = false;
        
        std::shared_ptr<rawfield> myfield;
        
    public:
        
        opfieldorder(std::shared_ptr<rawfield> fieldin) { myfield = fieldin; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
