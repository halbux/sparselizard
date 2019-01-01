// sparselizard - Copyright (C) 2017- A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
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
        shared_ptr<rawmat> rawmatptr = NULL;
        
        void errorifpointerisnull(void);
        	
	public:
                	
        mat(void) {};
        mat(shared_ptr<rawmat> inputrawmat) { rawmatptr = inputrawmat; };
     
    	// Create a matrix with structure based on a formulation and with initial values:
        mat(formulation myformulation, intdensematrix rowadresses, intdensematrix coladresses, densematrix vals);
     
        int countrows(void);
        int countcolumns(void);
        
        int countnnz(void);
        
        // Remove the rows and columns associated to Dirichlet constraints:
        void removeconstraints(void);
        
        void reuselu(void);
        
        shared_ptr<rawmat> getpointer(void) { return rawmatptr; };
        
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
