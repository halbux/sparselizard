// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <iostream>
#include <string>
#include <vector>
#include "expression.h"

class expression;

class integration
{
    private:
        
        std::vector<expression> myexpression = {};
        // Empty if undefined:
        std::vector<expression> mymeshdeform = {};
                
        int myphysicalregion;
        int myblocknumber;
        int myintegrationorderdelta;
        
        int mynumcoefharms = -1;

    public:
        
        // Barycenter eval mode for a rhs term:
        bool isbarycentereval = false;
        
        integration(void) {};
        
        integration(int physreg, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
        integration(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
        // For the multiharmonic computation an extra integer is required 
        // to know with how many harmonics the coef multiplying the tf 
        // and/or dof should be approximated. Set 'numcoefharms' to a 
        // negative integer and it will be as if you were calling 
        // the non multiharmonic functions above.
        integration(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
        integration(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);

        expression getexpression(void);
        bool ismeshdeformdefined(void) { return (mymeshdeform.size() == 1); };
        expression getmeshdeform(void);
        
        int getphysicalregion(void) { return myphysicalregion; };
        int getintegrationorderdelta(void) { return myintegrationorderdelta; };
        int getblocknumber(void) { return myblocknumber; };
        
        bool isfftrequested(void) { return (mynumcoefharms > 0); };
        int getnumberofcoefharms(void) { return mynumcoefharms; };
        
        void print(void);
};

#endif
