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

        int formfunctioncomponent = 0;

        // Cannnot be a 'x', 'y', 'z' or 'one' field:
        std::shared_ptr<rawfield> myfield;

    public:

        opfieldnosync(std::shared_ptr<rawfield> fieldin);

        void selectformfunctioncomponent(int comp) { formfunctioncomponent = comp; };

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);

};

#endif
