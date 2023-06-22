// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPSUM_H
#define OPSUM_H

#include "operation.h"

class opsum: public operation
{

    private:
        
        bool reuse = false;
        std::vector<std::shared_ptr<operation>> sumterms = {};
    
    public:
        
        opsum(void) {};
        opsum(std::vector<std::shared_ptr<operation>> input) { sumterms = input; };
        
        void addterm(std::shared_ptr<operation> term) { sumterms.push_back(term); }
        void subtractterm(std::shared_ptr<operation> term);
        void removeterm(int whichterm) { sumterms.erase(sumterms.begin() + whichterm); };
        
        int count(void) { return sumterms.size(); };
        bool issum(void) { return true; };

        std::vector<std::shared_ptr<operation>> getarguments(void) {return sumterms;};
        std::shared_ptr<operation> getargument(int argnum) { return sumterms[argnum]; };
        void replaceargument(int argnum, std::shared_ptr<operation> newarg) { sumterms[argnum] = newarg; };

        std::vector<std::vector<densemat>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densemat multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<operation> expand(void);
        // Group all sumterms and the sum terms they 
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
