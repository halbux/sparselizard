#ifndef OPPOWER_H
#define OPPOWER_H

#include "operation.h"

class oppower: public operation
{

	private:
        
        bool reuse = false;
        std::shared_ptr<operation> mybase;
        std::shared_ptr<operation> myexponent;
	
	public:
        
        oppower(std::shared_ptr<operation> base, std::shared_ptr<operation> exponent) { mybase = base; myexponent = exponent; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {mybase, myexponent}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
