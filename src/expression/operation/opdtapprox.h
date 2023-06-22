// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPDTAPPROX_H
#define OPDTAPPROX_H

#include "operation.h"

class opdtapprox: public operation
{

    private:

        int mydtorder;

        std::shared_ptr<operation> myarg;

        // Approximated dtx, dtdtx:
        double mydtx, mydtdtx;
        
        std::vector<double> mydtbkps = {};

    public:

        opdtapprox(int dtorder, std::shared_ptr<operation> arg, double initdtx, double initdtdtx);

        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {myarg}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void nextimpliciteuler(double tinit, double dt);
        void nextgenalpha(double beta, double gamma, double alphaf, double alpham, double tinit, double dt);

        void approvetimestep(void);

        void print(void);

};

#endif
