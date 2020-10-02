// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPHARMONIC_H
#define OPHARMONIC_H

#include "operation.h"

class opharmonic: public operation
{

    private:
        
        bool reuse = false;
        std::shared_ptr<operation> myarg;
        int myharmnum;
        int mynumfftharms;
        
    public:
        
        opharmonic(int harmnum, std::shared_ptr<operation> arg, int numfftharms) { myharmnum = harmnum; myarg = arg; mynumfftharms = numfftharms; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        bool isharmonicone(std::vector<int> disjregs) { return true; };

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {myarg}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
