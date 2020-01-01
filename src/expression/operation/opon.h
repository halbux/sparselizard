// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPON_H
#define OPON_H

#include "operation.h"
#include "expression.h"

class opon: public operation
{

    private:
        
        bool reuse = false;
        
        int myphysreg;
        bool myerrorifnotfound;
        // No shift if empty:
        std::vector<expression> mycoordshift = {};
        std::shared_ptr<operation> myarg;
        
    public:
        
        opon(int physreg, expression* coordshift, std::shared_ptr<operation> arg, bool errorifnotfound);
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void);
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };

        void print(void);

};

#endif
