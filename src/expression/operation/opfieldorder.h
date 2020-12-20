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
        
        double myalpha = -1.0;
        double mythreshold = 0.0;
        
        std::vector<std::shared_ptr<rawfield>> myfields = {};
        
    public:
        
        opfieldorder(std::vector<std::shared_ptr<rawfield>> fieldsin, double alpha = -1.0, double absthres = 0.0);
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
