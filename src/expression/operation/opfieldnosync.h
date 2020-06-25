// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPFIELDNOSYNC_H
#define OPFIELDNOSYNC_H

#include "operation.h"

class opfieldnosync: public operation
{

    private:

        int myformfunctioncomponent = -1;
        
        std::vector<std::shared_ptr<opfieldnosync>> mycomponents = {NULL};

        // Cannnot be a 'x', 'y', 'z' or 'one' field:
        std::shared_ptr<rawfield> myfield;

    public:

        opfieldnosync(int formfunctioncomponent, std::shared_ptr<rawfield> fieldin);
        // Provide the operations for all components after constructor:
        void setothercomponents(std::vector<std::shared_ptr<opfieldnosync>> allcomps);

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        
        void print(void);

};

#endif
