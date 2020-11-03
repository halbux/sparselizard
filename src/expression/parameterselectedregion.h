// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PARAMETERSELECTEDREGION_H
#define PARAMETERSELECTEDREGION_H

#include <iostream>
#include "expression.h"
#include "rawparameter.h"

class rawparameter;
class expression;

class parameterselectedregion
{

    private:
        
        std::shared_ptr<rawparameter> myparam;
        int myphysreg;
    
    public:
        
        parameterselectedregion(std::shared_ptr<rawparameter> param, int physreg) { myparam = param; myphysreg = physreg; };

        void operator=(expression input);
};

#endif
