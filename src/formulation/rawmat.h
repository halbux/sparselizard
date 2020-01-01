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
        
        long long int nnz = -1;
        
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
        
        std::shared_ptr<dofmanager> mydofmanager = NULL;
            
    public:
                    
        rawmat(std::shared_ptr<dofmanager> dofmngr) { mydofmanager = dofmngr; };
        rawmat(std::shared_ptr<dofmanager> dofmngr, Mat input) { mydofmanager = dofmngr; mymat = input; };

        ~rawmat(void);
     
        long long int countrows(void);
        long long int countcolumns(void);
        
        long long int countnnz(void) { return nnz; };
        
        // Set all row and/or column indices requested to -1 (-1 adresses are ignored at assembly):
        void zeroentries(intdensematrix entriestozero, bool zerorows, bool zerocolumns);

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
        // Create the petsc matrix.
        void process(void);
        // Remove the last added fragment:
        void removelastfragment(void);
        // Clear all the fragments:
        void clearfragments(void);
        
        void print(void);

        // Extract a new initialised rawmat that has all accumulated data. 
        std::shared_ptr<rawmat> extractaccumulated(void);
        
        std::shared_ptr<dofmanager> getdofmanager(void) { return mydofmanager; };
        Mat getpetsc(void) { return mymat; };
        
        KSP* getksp(void) { return &myksp; };

};

#endif
