// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object manages the data of a field, i.e. the coefficients in the
// finite element discretisation of the field: requested entries that
// are non existing are automatically added and filled with zeros.
// The field in this object can only contain a single component and 
// harmonic. The x, y and z coordinate fields are not supported.

#ifndef COEFMANAGER_H
#define COEFMANAGER_H

#include <string>
#include "disjointregions.h"
#include "hierarchicalformfunction.h"
#include <memory>
#include "selector.h"

class coefmanager
{

    private:

        disjointregions mydisjointregions;

        std::string myfieldtypename;

        // 'coefs[disjreg][formfunc][elem]' gives the coefficient for  
        //
        // - vertex, edge, face or volume type-disjoint region 'disjreg'
        // - the 'formfunc'th vertex, edge, face or volume form function
        // - element index 'elem' in the disjoint region
        //
        std::vector<std::vector<std::vector<double>>> coefs;

    public:

        coefmanager() {};
        coefmanager(std::string fieldtypename, disjointregions* drs);
        
        bool isdefined(int disjreg, int formfunctionindex);
        int countformfunctions(int disjreg);

        // Update the number of form functions considered in every disjoint region. 
        // To be called every time the interpolation order of the field changes.
        void fitinterpolationorder(int disjreg, int interpolationorder);

        double getcoef(int disjreg, int formfunctionindex, int elementindexindisjointregion);
        void setcoef(int disjreg, int formfunctionindex, int elementindexindisjointregion, double val);
        
        void print(bool databoundsonly);
        
};

#endif
