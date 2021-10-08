// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef H1POINT_H
#define H1POINT_H

#include "hierarchicalformfunction.h"

class h1point: public hierarchicalformfunction
{
    private:
    
        int targetdim = -1;

    public:
    
        h1point(int td) { targetdim = td; }; // -1 for h1

        // Get the number of form functions of order <= 'order':
        int count(int order);
        // Get the number of form functions of order <= 'order' that are associated to the num th
        // - node   in case dim = 0 
        // - edge   in case dim = 1
        // - face   in case dim = 2 
        // - volume in case dim = 3 
        int count(int order, int dim, int num);

        // Get the number of components in the form function.
        int countcomponents(void) { return 1; };

        // 'evalat' takes an integer giving the highest order up to which to output the form function polynomials.
        hierarchicalformfunctioncontainer evalat(int maxorder);

        // If 'isorientationdependent' is false then the assembly can
        // be carried out without taking care of the element orientation.
        // This provides an assembly speedup.
        bool isorientationdependent(int order) { return false; };
};

#endif
