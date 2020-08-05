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

        std::shared_ptr<operation> myarg;

        // Evaluation will be performed on this mesh state:
        std::shared_ptr<rawmesh> myrawmesh = NULL;
        std::shared_ptr<ptracker> myptracker = NULL;

    public:

        opnosync(std::shared_ptr<operation> arg, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt);

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        
        void print(void);

};

#endif
