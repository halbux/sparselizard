// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef OPCONSTANT_H
#define OPCONSTANT_H

#include "operation.h"

class opconstant: public operation
{

	private:
        
        bool reuse = false;
        double constantvalue;
	
	public:
        
        opconstant(double val) { constantvalue = val; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        bool isconstant(void) { return true; };
        double getvalue(void) { return constantvalue; };
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
