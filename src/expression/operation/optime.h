// sparselizard - Copyright (C) 2017- A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This object treats the time variable 't'.

#ifndef OPTIME_H
#define OPTIME_H

#include "operation.h"

class optime: public operation
{

	private:
        
        bool reuse = false;
        
	public:
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
