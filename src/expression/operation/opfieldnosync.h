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

        int mycomp = -1;
        
        std::vector<std::weak_ptr<opfieldnosync>> mycomponents = {};

        std::shared_ptr<rawfield> myfield = NULL;

    public:

        opfieldnosync(int formfunctioncomponent, std::shared_ptr<rawfield> fieldin);
        // Provide the operations for all components after constructor:
        void setcomponents(std::vector<std::shared_ptr<opfieldnosync>> allcomps);

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        
        void print(void);

};

#endif
