// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPPORT_H
#define OPPORT_H

#include "operation.h"
#include "rawport.h"

class opport: public operation
{

    private:
        
        int timederivativeorder = 0;

        std::shared_ptr<rawport> myport;
    
    public:
        
        opport(std::shared_ptr<rawport> portin);

        void increasetimederivativeorder(int amount);

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        bool isport(void) { return true; };

        int gettimederivative(void) { return timederivativeorder; };

        std::shared_ptr<rawport> getportpointer(void);

        bool isharmonicone(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);

        void print(void);

};

#endif
