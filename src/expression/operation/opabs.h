#ifndef OPABS_H
#define OPABS_H

#include "operation.h"

class opabs: public operation
{

	private:
        
        bool reuse = false;
        std::shared_ptr<operation> myarg;
        
	public:
        
        opabs(std::shared_ptr<operation> arg) { myarg = arg; };
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {myarg}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        
        void print(void);

};

#endif
