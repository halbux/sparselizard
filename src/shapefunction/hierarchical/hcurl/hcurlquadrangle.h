// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// The following form functions have been derived in the PhD thesis
// "High Order Finite Element Methods for Electromagnetic Field Computation" 
// of Dr. Sabine Zaglmayr at JOHANNES KEPLER UNIVERSITAET LINZ, AUSTRIA

#ifndef HCURLQUADRANGLE_H
#define HCURLQUADRANGLE_H

#include "hierarchicalformfunction.h"

class hcurlquadrangle: public hierarchicalformfunction
{
    private:

    public:
    
        // Get the number of form functions of order <= 'order':
        int count(int order);
        // Get the number of form functions of order <= 'order' that are associated to the num th
        // - node   in case dim = 0 
        // - edge   in case dim = 1
        // - face   in case dim = 2 
        // - volume in case dim = 3 
        int count(int order, int dim, int num);

        // Get the number of components in the form function.
        int countcomponents(void) { return 3; };

        // 'evalat' takes an integer giving the highest order up to which to output the form function polynomials.
        hierarchicalformfunctioncontainer evalat(int maxorder);

        // Return a vector whose index i is true if the ith form function
        // associated to highest dimension elements is of gradient type.
        std::vector<bool> isgradienttype(int order);
};

#endif
