// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPCUSTOM_H
#define OPCUSTOM_H

#include "operation.h"

class opcustom: public operation
{

    private:
        
        // Custom operations are always reused
        
        int myoutindex = -1;
        
        std::vector<std::shared_ptr<operation>> myargs = {};
        std::vector<field> myfields = {}; // for advanced custom function
        
        std::vector<std::weak_ptr<opcustom>> myfamily = {};
        
        std::vector<densemat> (*myfunction)(std::vector<densemat>) = NULL;
        std::vector<densemat> (*myadvancedfunction)(std::vector<densemat>, std::vector<field>, elementselector&, std::vector<double>&, expression*) = NULL;
        
    public:
        
        opcustom(int outindex, std::vector<densemat> fct(std::vector<densemat>), std::vector<std::shared_ptr<operation>> args);
        opcustom(int outindex, std::vector<densemat> fct(std::vector<densemat>, std::vector<field>, elementselector&, std::vector<double>&, expression*), std::vector<std::shared_ptr<operation>> args, std::vector<field> infields);
        
        // Provide all related operations:
        void setfamily(std::vector<std::weak_ptr<opcustom>> ops) { myfamily = ops; };
        
        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return myargs; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void print(void);

};

#endif
