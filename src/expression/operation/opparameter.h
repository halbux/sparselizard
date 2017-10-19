// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef OPPARAMETER_H
#define OPPARAMETER_H

#include "operation.h"

class parameter;

class opparameter: public operation
{

	private:
        
        bool reuse = false;
        
        int myrow;
        int mycolumn;
        
        parameter* myparameter;
	
	public:
        
        opparameter(parameter* input, int row, int col) { myparameter = input; myrow = row; mycolumn = col; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        
        parameter* getparameterpointer(void) { return myparameter; };
        
        bool isharmonicone(std::vector<int> disjregs);
        
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        bool isvalueorientationdependent(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
