// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef ONECONSTANT_H
#define ONECONSTANT_H

#include "hierarchicalformfunction.h"

class oneconstant: public hierarchicalformfunction
{
    private:
    
        int targetdim = -1;
        int elementtypenumber = -1;
        int elementdimension = -1;

    public:
    
        oneconstant(int td, int et);

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
