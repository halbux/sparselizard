// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object is a wrapper of the actual vec object 'rawvec' pointed
// by 'rawvecptr' and to which all functions are redirected.

#ifndef VEC_H
#define VEC_H

#include "formulation.h"
#include <iostream>
#include <string>
#include <numeric>
#include "petsc.h"
#include "field.h"
#include "vectorfieldselect.h"
#include "memory.h"
#include "rawvec.h"
#include "petscvec.h"

class field;
class rawvec;
class formulation;

class vec
{
    private:

        // The actual vector:
        shared_ptr<rawvec> rawvecptr = NULL;
        
        void errorifpointerisnull(void);
    
    public:
            
        vec(void) {};
        vec(formulation formul);
        vec(shared_ptr<rawvec> inputrawvec) { rawvecptr = inputrawvec; };
     
        int size(void);
        
        // Remove the entries associated to Dirichlet constraints:
        void removeconstraints(void);
        
        // Update the Dirichlet constraints:
        void updateconstraints(void);
        
        // Negative addresses are ignored. 'op' can be 'add' or 'set'. 
        void setvalues(intdensematrix addresses, densematrix valsmat, std::string op = "set");
        densematrix getvalues(intdensematrix addresses);
        // Set and get value at a single index:
        void setvalue(int address, double value, std::string op = "set");
        double getvalue(int address);
        
        vectorfieldselect operator|(field selectedfield);
        
        shared_ptr<rawvec> getpointer(void) { return rawvecptr; };
        
        // Transfer data from a field to this vector.
        // Only the data corresponding to the physical region is transferred.
        // 'op' can be 'add' or 'set'.
        void setdata(int physreg, field myfield, std::string op = "set");
        
        Vec getpetsc(void);
        
        // Write and load raw vec data:
        void write(std::string filename);
        void load(std::string filename);
        
        void print(void);
        
        vec copy(void);
        
        double norm(std::string type = "2");
        double sum(void);
        
        vec operator+(void);
        vec operator-(void);
        vec operator*(double input);
        vec operator+(vec input);
        vec operator-(vec input);
        
};

vec operator*(double, vec);

#endif
