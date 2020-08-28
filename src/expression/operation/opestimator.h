// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPESTIMATOR_H
#define OPESTIMATOR_H

#include "operation.h"
#include "opfield.h"

class opfield;

class opestimator: public operation
{

    private:
        
        bool reuse = false;
        // Estimator type:
        std::string mytype = "";
        std::shared_ptr<operation> myarg;
        
        // This contains the estimated value:
        std::shared_ptr<opfield> myvalue = NULL;
        
        long long int mystatenumber = 0;
        
    public:
        
        opestimator(std::string estimatortype, std::shared_ptr<operation> arg);
        
        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::vector<std::shared_ptr<operation>> getarguments(void) { return {myarg}; };
        std::shared_ptr<operation> simplify(std::vector<int> disjregs);
        
        bool isvalueorientationdependent(std::vector<int> disjregs) { return false; };
        
        std::shared_ptr<operation> copy(void);
        
        void reuseit(bool istobereused) { reuse = istobereused; };

        void print(void);
        
        // Update 'myvalue' with the estimate:
        void estimatezienkiewiczzhu(void);

};

#endif
