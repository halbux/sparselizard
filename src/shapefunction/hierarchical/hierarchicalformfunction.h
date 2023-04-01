// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef HIERARCHICALFORMFUNCTION_H
#define HIERARCHICALFORMFUNCTION_H

#include "math.h"
#include <algorithm>
#include "orientation.h"
#include "element.h"
#include "legendre.h"
#include "polynomial.h"
#include <vector>
#include <string>
#include "hierarchicalformfunctioncontainer.h"

class hierarchicalformfunction
{
    private:
    
        // List all form function type names (type numbers are defined by the name ordering):
        std::vector<std::string> mytypenames = {"h1","hcurl","one0","one1","one2","one3","h1d0","h1d1","h1d2","h1d3"};

    public:

        // Translate a form function type number to a type name and vice versa:
        int gettypenumber(std::string fftypename);
        std::string gettypename(int fftypenumber);
        
        // Get the minimum order allowed:
        int getminorder(std::string fftypename);

        // Get the number of form functions of order <= 'order':
        virtual int count(int order) { throw std::runtime_error(""); }; // fix return warning
        // Get the number of form functions of order <= 'order' that are associated to the num th
        // - node   in case dim = 0 
        // - edge   in case dim = 1
        // - face   in case dim = 2 
        // - volume in case dim = 3 
        virtual int count(int order, int dim, int num) { throw std::runtime_error(""); }; // fix return warning

        // Get the number of components in the form function.
        virtual int countcomponents(void) { throw std::runtime_error(""); }; // fix return warning

        // 'evalat' takes an integer giving the highest order up to which to output the form function polynomials.
        virtual hierarchicalformfunctioncontainer evalat(int maxorder) { throw std::runtime_error(""); }; // fix return warning

        // If 'isorientationdependent' is false then the assembly can
        // be carried out without taking care of the element orientation.
        // This provides an assembly speedup. By default it is not used.
        virtual bool isorientationdependent(int order) { return true; };

        // Return a vector whose index i is true if the ith form function
        // associated to highest dimension elements is of gradient type.
        virtual std::vector<bool> isgradienttype(int order) { return std::vector<bool>(0); };
};

#endif
