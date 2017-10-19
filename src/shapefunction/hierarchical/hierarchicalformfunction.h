// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef HIERARCHICALFORMFUNCTION_H
#define HIERARCHICALFORMFUNCTION_H

#include "math.h"
#include <algorithm>
#include "orientation.h"
#include "element.h"
#include "legendre.h"
#include "polynomial.h"
#include <vector>
#include "hierarchicalformfunctioncontainer.h"

class hierarchicalformfunction
{
	private:
	
	public:
	
        // Get the number of form functions of order <= 'order':
        virtual int count(int order) {};
        // Get the number of form functions of order <= 'order' that are associated to the num th
        // - node   in case dim = 0 
        // - edge   in case dim = 1
        // - face   in case dim = 2 
        // - volume in case dim = 3 
        virtual int count(int order, int dim, int num) {};
        
        // Get the number of components in the form function.
        virtual int countcomponents(void) {};

        // 'evalat' takes a vector with the coordinates of the evaluation 
        // points in the format [ki1 eta1 phi1 ki2 eta2 phi2 ki3 ...] as well 
        // as an integer giving the highest order up to which to output the form 
        // function polynomials.
        virtual hierarchicalformfunctioncontainer evalat(int maxorder, std::vector<double> evaluationpoints) {};

        // If 'isorientationdependent' is false then the assembly can
        // be carried out without taking care of the element orientation.
        // This provides an assembly speedup. By default it is not used.
        virtual bool isorientationdependent(int order) { return true; };
        
};

#endif
