// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef ONCONTEXT_H
#define ONCONTEXT_H

#include <iostream>
#include <vector>
#include <string>
#include "expression.h"

class oncontext
{

    private:

        bool myisdefined = false;    
    
        int myphysreg = -1;
        bool myerrorifnotfound = true;
        // Empty if not shifted:
        std::vector<expression> mycoordshift = {};
    
    public:

        oncontext(void) {};
        oncontext(int physreg, expression* coordshift, bool errorifnotfound);            

        bool isdefined(void);

        int getphysicalregion(void);
        
        bool isshifted(void);
        expression* getshift(void);
        bool iserrorifnotfound(void);
        
        // Compare two oncontexes:
        bool isequal(oncontext* tocompare);
            
};

#endif

