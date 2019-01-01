// sparselizard - Copyright (C) 2017- A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef HCURLPOINT_H
#define HCURLPOINT_H

#include "hierarchicalformfunction.h"

class hcurlpoint: public hierarchicalformfunction
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

        // 'evalat' takes a vector with the coordinates of the evaluation 
        // points in the format [ki1 eta1 phi1 ki2 eta2 phi2 ki3 ...] as well 
        // as an integer giving the highest order up to which to output the form 
        // function polynomials.
        hierarchicalformfunctioncontainer evalat(int maxorder, std::vector<double> evaluationpoints);
};

#endif
