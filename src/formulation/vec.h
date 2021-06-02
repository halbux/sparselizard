// sparselizard - Copyright (C) see copyright file.
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
#include "port.h"

class field;
class rawvec;
class formulation;

class vec
{
    private:

        // The actual vector:
        std::shared_ptr<rawvec> rawvecptr = NULL;
        
        void errorifpointerisnull(void);
    
    public:
            
        vec(void) {};
        vec(formulation formul);
        vec(std::shared_ptr<rawvec> inputrawvec) { rawvecptr = inputrawvec; };
     
        // Create a pre-filled vector:
        vec(int vecsize, intdensematrix addresses, densematrix vals);
     
        int size(void);
        
        void permute(intdensematrix rowpermute, bool invertit = false);
        
        // Update the Dirichlet constraints:
        void updateconstraints(void);
        
        // Negative addresses are ignored. 'op' can be 'add' or 'set'. 
        void setvalues(intdensematrix addresses, densematrix valsmat, std::string op = "set");
        void setallvalues(densematrix valsmat, std::string op = "set");
        densematrix getvalues(intdensematrix addresses);
        densematrix getallvalues(void);
        // Set and get value at a single index:
        void setvalue(int address, double value, std::string op = "set");
        double getvalue(int address);
        // Set and get value of a port:
        void setvalue(port prt, double value, std::string op = "set");
        double getvalue(port prt);
        
        vectorfieldselect operator|(field selectedfield);
        
        std::shared_ptr<rawvec> getpointer(void) { return rawvecptr; };
        
        // Transfer data from a field to this vector.
        // Only the data corresponding to the physical region is transferred.
        // 'op' can be 'add' or 'set'.
        void setdata(int physreg, field myfield, std::string op = "set");
        // Transfer data from all fields and ports to this vector.
        void setdata(void);
        
        // Allow/forbid automatic updating of the vec value during hp-adaptivity:
        void automaticupdate(bool updateit);
        void noautomaticupdate(void);
        
        Vec getpetsc(void);
        
        // Write and load raw vec data:
        void write(std::string filename);
        void load(std::string filename);
        
        void print(void);
        
        vec copy(void);
        
        vec extract(intdensematrix addresses);
        
        double norm(std::string type = "2");
        double sum(void);
        
        vec operator+(void);
        vec operator-(void);
        vec operator*(double input);
        vec operator/(double input);
        vec operator+(vec input);
        vec operator-(vec input);
        
};

vec operator*(double, vec);

#endif
