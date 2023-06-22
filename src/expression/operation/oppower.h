// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPPOWER_H
#define OPPOWER_H

#include "operation.h"

class oppower: public operation
{

    private:
        
        bool reuse = false;
        std::shared_ptr<operation> mybase;
        std::shared_ptr<operation> myexponent;
    
    public:
        
        oppower(std::shared_ptr<operation> base, std::shared_ptr<operation> exponent) { mybase = base; myexponent = exponent; };
        
        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {mybase, myexponent}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);

        void print(void);

};

#endif
