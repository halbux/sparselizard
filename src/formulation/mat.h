// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object is a wrapper of the actual mat object 'rawmat' pointed
// by 'rawmatptr' and to which all functions are redirected.

#ifndef MAT_H
#define MAT_H

#include <iostream>
#include "petsc.h"
#include <vector>
#include "memory.h"
#include "rawmat.h"
#include "petscvec.h" 
#include "petscmat.h"
#include "vec.h"
#include "rawvec.h"
#include "formulation.h"
#include "intdensematrix.h"
#include "densematrix.h"

class vec;
class rawmat;
class formulation;

class mat
{
    private:
        
        // The actual matrix:
        std::shared_ptr<rawmat> rawmatptr = NULL;
        
        void errorifpointerisnull(void);
        void errorifinvalidated(void);
            
    public:
                    
        mat(void) {};
        mat(std::shared_ptr<rawmat> inputrawmat) { rawmatptr = inputrawmat; };
        
        // Create a pre-filled matrix:
        mat(long long int matsize, intdensematrix rowadresses, intdensematrix coladresses, densematrix vals);
        // Create a matrix with structure based on a formulation and with initial values:
        mat(formulation myformulation, intdensematrix rowadresses, intdensematrix coladresses, densematrix vals);
     
        bool isdefined(void) { return (rawmatptr != NULL); };
         
        long long int countrows(void);
        long long int countcolumns(void);
        
        long long int countnnz(void);
        
        // Permute the rows and columns:
        void permute(intdensematrix rowpermute, intdensematrix colpermute);
        
        // Remove the rows and columns associated to Dirichlet constraints:
        void removeconstraints(void);
        
        void reusefactorization(void);
        
        std::shared_ptr<rawmat> getpointer(void);
        
        Mat getpetsc(void);
        
        void print(void);
        
        mat copy(void);
        
        
        mat operator+(void);
        mat operator-(void);
        mat operator*(double input);
        mat operator*(mat input);
        mat operator+(mat input);
        mat operator-(mat input);
        vec operator*(vec input);
        
};

mat operator*(double, mat);

#endif
