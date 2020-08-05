// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPNOSYNC_H
#define OPNOSYNC_H

#include "operation.h"

class opnosync: public operation
{

    private:

        int mycomp = -1;
        
        std::vector<std::weak_ptr<opnosync>> mycomponents = {};

        std::shared_ptr<rawfield> myfield = NULL;

    public:

        opnosync(int formfunctioncomponent, std::shared_ptr<rawfield> fieldin);
        // Provide the operations for all components after constructor:
        void setcomponents(std::vector<std::shared_ptr<opnosync>> allcomps);

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        
        void print(void);

};

#endif
