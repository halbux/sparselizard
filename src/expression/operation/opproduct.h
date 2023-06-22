// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPPRODUCT_H
#define OPPRODUCT_H

#include "operation.h"
#include "harmonic.h"
#include "fourier.h"

class opproduct: public operation
{

    private:
        
        bool reuse = false;
        std::vector<std::shared_ptr<operation>> productterms = {};
    
    public:
        
        opproduct(void) {};
        opproduct(std::vector<std::shared_ptr<operation>> input) { productterms = input; };
                
        void multiplybyterm(std::shared_ptr<operation> term) { productterms.push_back(term); };
        
        void removeterm(int whichterm) { productterms.erase(productterms.begin() + whichterm); };
        
        int count(void) { return productterms.size(); };
        bool isproduct(void) { return true; };

        std::vector<std::shared_ptr<operation>> getarguments(void) {return productterms;};
        std::shared_ptr<operation> getargument(int argnum) { return productterms[argnum]; };
        void replaceargument(int argnum, std::shared_ptr<operation> newarg) { productterms[argnum] = newarg; };
        
        // Multiharmonic products are not allowed unless one of 
        // the two product arguments has only the cos0 harmonic.
        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> expand(void);
        // Group all productterms and the product terms they 
        // include recursively together in this object.
        void group(void);
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };
        bool isreused(void) { return reuse; };
        
        std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);

        void print(void);

};

#endif
