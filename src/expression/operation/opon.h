// sparselizard - Copyright (C) 2017- A. Halbach
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
        // No shift if empty:
        std::vector<expression> mycoordshift = {};
        std::shared_ptr<operation> myarg;
        
	public:
        
        opon(int physreg, expression* coordshift, std::shared_ptr<operation> arg);
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void);
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
		/// REMOVE + CHECK ERROR OK IF USED std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);

        void print(void);

};

#endif
