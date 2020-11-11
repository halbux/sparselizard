// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This function computes (at the given Gauss points) the Jacobian of the 
// variable change bringing back the mesh element to the reference element.
// If requested the calculations are carried out on the mesh deformed by
// the provided vector expression 'meshdeform'. The expression must be 
// constant in time (harmonic 1).
//
// The Jacobian matrix is defined as follows:
//
//           dx/dki   dy/dki   dz/dki
//           dx/deta  dy/deta  dz/deta
//           dx/dphi  dy/dphi  dz/dphi
//
// Its inverse therefore equals:
//
//           dki/dx   deta/dx  dphi/dx
//           dki/dy   deta/dy  dphi/dy
//           dki/dz   deta/dz  dphi/dz
//

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <iostream>
#include <vector>
#include <math.h>
#include <memory>
#include "field.h"
#include "rawfield.h"
#include "expression.h"
#include "universe.h"
#include "densematrix.h"
#include "elementselector.h"
#include "disjointregions.h"

class field;
class expression;
class elementselector;

class jacobian
{

    private:
        
        densematrix detjac, xcoord;
        std::vector<densematrix> jac = std::vector<densematrix>(3*3);
        std::vector<densematrix> invjac = {};

    public:
        
        jacobian(void) {};
        
        // The Jacobian is computed on the mesh deformed by vector expression 
        // 'meshdeform'. The expression must be constant in time (harmonic 1).
        jacobian(elementselector& elemselect, std::vector<double> evaluationcoordinates, expression* meshdeform);

        // Return a jacobian object holding an element subset of this jacobian:
        jacobian extractsubset(std::vector<int>& selectedelementindexes);
        
        // The detjac and jac terms are computed in the constructor.
        densematrix getdetjac(void);    
        densematrix getjac(int row, int column);
        // The invjac terms are only computed when 'getinvjac' is called.
        densematrix getinvjac(int row, int column);    

};

#endif
