// sparselizard - Copyright (C) 2017- A. Halbach, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef PARAMETERSELECTEDREGION_H
#define PARAMETERSELECTEDREGION_H

#include <iostream>
#include "expression.h"
#include "parameter.h"

class parameter;
class expression;

class parameterselectedregion
{

	private:
        
        parameter* myparam;
        int myphysreg;
	
	public:
        
        parameterselectedregion(parameter* param, int physreg) { myparam = param; myphysreg = physreg; };

        void operator=(expression input);
};

#endif
