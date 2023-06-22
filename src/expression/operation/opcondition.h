// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPCONDITION_H
#define OPCONDITION_H

#include "operation.h"

class opcondition: public operation
{

    private:
        
        bool reuse = false;
        std::shared_ptr<operation> mycond;
        std::shared_ptr<operation> mytrue;
        std::shared_ptr<operation> myfalse;
        
    public:
        
        opcondition(std::shared_ptr<operation> condarg, std::shared_ptr<operation> truearg, std::shared_ptr<operation> falsearg);
        
        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {mycond, mytrue, myfalse}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);

        void print(void);

};

#endif
