// sparselizard - Copyright (C) 2017- A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
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
	
		shared_ptr<dofmanager> mydofmanager;
        
        // All elementary coef*dof*tf terms to assemble and add together:
        std::vector<shared_ptr<operation>> mydofs = {};
        std::vector<shared_ptr<operation>> mytfs = {};
        std::vector<shared_ptr<operation>> mycoeffs = {};
        
        // The dof and tf field for all terms above. A NULL dof means rhs contribution:
        shared_ptr<rawfield> doffield = NULL;
        shared_ptr<rawfield> tffield = NULL;

        int integrationphysreg = -1;
        int dofphysreg = -1;
        int tfphysreg = -1;
        
        int integrationorderdelta = 0;
        // Number of time evaluations for the FFT of the coef. Negative means no FFT.
        int numfftcoeffs = -1;
        
        // The contribution is computed on the mesh deformed by (if any):
        std::vector<expression> mymeshdeformation = {};

        // The following vectors stores pointers to all the fragments that have been generated.
        std::vector<intdensematrix> fragmentrowadresses = {};
        std::vector<intdensematrix> fragmentcoladresses = {};
        std::vector<densematrix> fragmentvalues = {};

	public:
	
		contribution(shared_ptr<dofmanager> dofmngr);
        
        void setdofs(std::vector<shared_ptr<operation>> dofs);
        void settfs(std::vector<shared_ptr<operation>> tfs);
        void setcoeffs(std::vector<shared_ptr<operation>> coeffs);
        void setdoffield(shared_ptr<rawfield> input);
        void settffield(shared_ptr<rawfield> input);
        void setmeshdeformation(expression meshdeform);
        void setintegrationphysicalregion(int physreg);
        void setdofphysicalregion(int physreg);
        void settfphysicalregion(int physreg);
        void setintegrationorderdelta(int integrorderdelta);
        void setnumfftcoeffs(int numcoeffs);
        
        // Generate the contribution and store it in the 
        // vec (for rhs contributions) or in the mat.
		void generate(shared_ptr<rawvec> myvec, shared_ptr<rawmat> mymat, bool computeconstraints = true);
		
// Add this to reduce the size of the system to solve:
//void removebubblemode(void);
																	
};

#endif
