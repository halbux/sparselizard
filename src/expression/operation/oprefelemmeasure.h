// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.

// This object gives the length of the reference line, surface of the reference 
// triangle/quadrangle and volume of the reference volume elements.


#ifndef OPREFELEMMEASURE_H
#define OPREFELEMMEASURE_H

#include "operation.h"

class oprefelemmeasure: public operation
{

	private:
        
	public:
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> copy(void);
        
        void print(void);

};

#endif
