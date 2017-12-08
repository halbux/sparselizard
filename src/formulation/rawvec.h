// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.

// This code calls the PETSc library. See https://www.mcs.anl.gov/petsc/ for more information.

#ifndef RAWVEC_H
#define RAWVEC_H

#include <iostream>
#include <string>
#include "dofmanager.h"
#include "intdensematrix.h"
#include "densematrix.h"
#include "vectorfieldselect.h"
#include "rawfield.h"
#include "memory.h"
#include "petsc.h"
#include "petscvec.h"
#include "mathop.h"

class vectorfieldselect;
class dofmanager;
class rawfield;

class rawvec
{
	private:

        Vec myvec;
        shared_ptr<dofmanager> mydofmanager = NULL;
	
	public:
        	
        rawvec(shared_ptr<dofmanager> dofmngr);
        rawvec(shared_ptr<dofmanager> dofmngr, Vec input) { mydofmanager = dofmngr; myvec = input; };
        
        ~rawvec(void);
        
        int size(void);
        
        // Remove the entries associated to Dirichlet constraints:
        void removeconstraints(void);
        
        // Update the indexes that correspond to constrained 
        // values of a rawfield on given disjoint regions. 
        void updateconstraints(shared_ptr<rawfield> constrainedfield, std::vector<int> disjregs);
     
        // Negative adresses are ignored. 'op' can be 'add' or 'set'. 
		void setvalues(intdensematrix adresses, densematrix valsmat, std::string op = "set");

		densematrix getvalues(intdensematrix adresses);
		densematrix getvalues(shared_ptr<rawfield> selectedfield, int disjointregionnumber, int formfunctionindex);
        
		void print(void);

        shared_ptr<dofmanager> getdofmanager(void) { return mydofmanager; };
        Vec getpetsc(void) { return myvec; };

		// Tranfer data between the 'inputvec' vector and this vector:
        void setdata(shared_ptr<rawvec> inputvec, int disjreg, shared_ptr<rawfield> inputfield);
        
};

#endif
