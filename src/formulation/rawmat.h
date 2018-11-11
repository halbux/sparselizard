// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.

// This code calls the PETSc library. See https://www.mcs.anl.gov/petsc/ for more information.

#ifndef RAWMAT_H
#define RAWMAT_H

#include <iostream>
#include <vector>
#include "dofmanager.h"
#include <cmath>
#include "myalgorithm.h"
#include "intdensematrix.h"
#include "densematrix.h"
#include "memory.h"
#include "petsc.h"
#include "petscmat.h"

class dofmanager;

class rawmat
{
	private:
        
        int nnz = -1;
        
        // Combining all accumulated fragments below gives the overall matrix.
        std::vector<intdensematrix> accumulatedrowindices = {};
        std::vector<intdensematrix> accumulatedcolindices = {};
        std::vector<densematrix> accumulatedvals = {};
        
        // Petsc does not deallocate the matrix data, we have to keep track of it.
        intdensematrix petscrows, petsccols;
        densematrix petscvals;

        Mat mymat = PETSC_NULL;
        
        // 'myksp' will store the LU decomposition if it is to be reused:
        KSP myksp = PETSC_NULL;
        bool lureuse = false;
        bool ludefined = false;
        
        shared_ptr<dofmanager> mydofmanager = NULL;
        	
	public:
                	
        rawmat(shared_ptr<dofmanager> dofmngr) { mydofmanager = dofmngr; };
        rawmat(shared_ptr<dofmanager> dofmngr, Mat input) { mydofmanager = dofmngr; mymat = input; };

        ~rawmat(void);
     
        int countrows(void);
        int countcolumns(void);
        
        int countnnz(void) { return nnz; };

        // Set the gauged indices to -1:
        void gauge(void);  
        
        // Remove the rows and columns associated to Dirichlet constraints:
        void removeconstraints(void);
        
        void reuselu(void) { lureuse = true; };
        bool islutobereused(void) { return lureuse; };
        bool isludefined(void) { return ludefined; };
        void isludefined(bool isdefined) { ludefined = isdefined; };
    
        // Add a fragment to the matrix.
        void accumulate(intdensematrix rowadresses, intdensematrix coladresses, densematrix vals);   
        // Create the petsc matrix. Choose to keep the fragments for reuse or not.
        void process(bool keepfragments = false);
        // Remove the last added fragment:
        void removelastfragment(void);
        
        void print(void);

    	// Extract a new initialised rawmat that has all accumulated data. 
    	shared_ptr<rawmat> extractaccumulated(void);
    	
        shared_ptr<dofmanager> getdofmanager(void) { return mydofmanager; };
        Mat getpetsc(void) { return mymat; };
        
        KSP* getksp(void) { return &myksp; };

};

#endif
