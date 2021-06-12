// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef CONTRIBUTION_H
#define CONTRIBUTION_H

#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include "operation.h"
#include "densematrix.h"
#include "intdensematrix.h"
#include "dofmanager.h"
#include "rawfield.h"
#include "universe.h"
#include "jacobian.h"
#include "gausspoints.h"
#include "elementselector.h"
#include "disjointregions.h"
#include "disjointregionselector.h"
#include "expression.h"
#include "harmonic.h"
#include "myfft.h"
#include "selector.h"
#include "math.h"
#include "rawvec.h"
#include "rawmat.h"
#include "wallclock.h"
#include "operation.h"

class rawvec;
class rawmat;
class operation;
class rawfield;

class contribution
{

    private:
    
        std::shared_ptr<dofmanager> mydofmanager;
        
        // All elementary coef*dof*tf terms to assemble and add together:
        std::vector<std::shared_ptr<operation>> mydofs = {};
        std::vector<std::shared_ptr<operation>> mytfs = {};
        std::vector<std::shared_ptr<operation>> mycoeffs = {};
        
        // The dof and tf field for all terms above. A NULL dof means rhs contribution:
        std::shared_ptr<rawfield> doffield = NULL;
        std::shared_ptr<rawfield> tffield = NULL;

        int integrationphysreg = -1;
        int dofphysreg = -1;
        int tfphysreg = -1;
        
        int integrationorderdelta = 0;
        // Number of time evaluations for the FFT of the coef. Negative means no FFT.
        int numfftcoeffs = -1;
        
        // Barycenter evaluation mode during rhs term assembly:
        bool isbarycentereval = false;
        
        // The contribution is computed on the mesh deformed by (if any):
        std::vector<expression> mymeshdeformation = {};

        // The following vectors stores pointers to all the fragments that have been generated.
        std::vector<intdensematrix> fragmentrowadresses = {};
        std::vector<intdensematrix> fragmentcoladresses = {};
        std::vector<densematrix> fragmentvalues = {};

    public:
    
        contribution(std::shared_ptr<dofmanager> dofmngr);
        
        void setdofs(std::vector<std::shared_ptr<operation>> dofs);
        void settfs(std::vector<std::shared_ptr<operation>> tfs);
        void setcoeffs(std::vector<std::shared_ptr<operation>> coeffs);
        void setdoffield(std::shared_ptr<rawfield> input);
        void settffield(std::shared_ptr<rawfield> input);
        void setmeshdeformation(expression meshdeform);
        void setintegrationphysicalregion(int physreg);
        void setdofphysicalregion(int physreg);
        void settfphysicalregion(int physreg);
        void setintegrationorderdelta(int integrorderdelta);
        void setnumfftcoeffs(int numcoeffs);
        void setbarycenterevalflag(void);
        
        // Generate the contribution and store it in the 
        // vec (for rhs contributions) or in the mat.
        void generate(std::shared_ptr<rawvec> myvec, std::shared_ptr<rawmat> mymat);
                                            
};

#endif
