// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

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
        
        // Combining all accumulated fragments below gives the overall matrix.
        std::vector<intdensematrix> accumulatedrowindices = {};
        std::vector<intdensematrix> accumulatedcolindices = {};
        std::vector<densematrix> accumulatedvals = {};
        

        long long int nnzA = -1, nnzD = -1;
        
        // The sparse matrix is stored in csr format. Dirichlet constraints are eliminated using A and D in Atotal = [A D; 0 1].
        // The rows and columns in A and D are renumbered consecutively from zero.
        intdensematrix Arows, Acols, Drows, Dcols;
        densematrix Avals, Dvals;
        // Ainds[i] is the index in Atotal of the ith index in A:
        intdensematrix Ainds, Dinds;
        
        Mat Amat = PETSC_NULL, Dmat = PETSC_NULL;
        
        // 'myksp' will store the factorization if it is to be reused:
        KSP myksp = PETSC_NULL;
        bool factorizationreuse = false;
        bool isitfactored = false;
        
        std::shared_ptr<dofmanager> mydofmanager = NULL;
        
        int mymeshnumber = 0;
            
    public:
                    
        rawmat(std::shared_ptr<dofmanager> dofmngr);
        rawmat(std::shared_ptr<dofmanager> dofmngr, Mat inA, Mat inD, intdensematrix inAinds, intdensematrix inDinds);

        ~rawmat(void);
     
        long long int countrows(void);
        long long int countcolumns(void);
        
        long long int countnnz(void) { return nnzA; };
        
        int getmeshnumber(void) { return mymeshnumber; };
        
        void reusefactorization(void) { factorizationreuse = true; };
        bool isfactorizationreuseallowed(void) { return factorizationreuse; };
        bool isfactored(void) { return isitfactored; };
        void isfactored(bool isfact) { isitfactored = isfact; };
    
        // Add a fragment to the matrix (empty fragments are ignored):
        void accumulate(intdensematrix rowadresses, intdensematrix coladresses, densematrix vals);   
        // Create the petsc matrices:
        void process(std::vector<bool>& isconstrained);
        // Clear all the fragments:
        void clearfragments(void);
        
        void print(void);

        // Extract a new initialized rawmat that has all accumulated data:
        std::shared_ptr<rawmat> extractaccumulated(void);
        
        intdensematrix getainds(void);
        intdensematrix getdinds(void);

        Mat getapetsc(void);
        Mat getdpetsc(void);
        
        std::shared_ptr<dofmanager> getdofmanager(void);
        
        KSP* getksp(void);

};

#endif
