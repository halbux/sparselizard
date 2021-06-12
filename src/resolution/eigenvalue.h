// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// Thanks to R. Haouari for the damped eigenvalue calculation code.

// This object uses the SLEPc library to get the eigenvalues and eigenvectors.
// More information on SLEPc can be found at http://slepc.upv.es

#ifndef EIGENVALUE_H
#define EIGENVALUE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include "vec.h"
#include "mat.h"
#include "rawvec.h"
#include "memory"

class eigenvalue
{
    private:
        
        mat myA;
        mat myB;
        std::vector<mat> mymats = {};
        
        // Real and imaginary part of the eigenvalues and eigenvectors:
        std::vector<double> eigenvaluereal = {};
        std::vector<double> eigenvalueimaginary = {};
        
        std::vector<vec> eigenvectorreal = {};
        std::vector<vec> eigenvectorimaginary = {};
        
    public:

        // Define a standard eigenvalue problem A*x = lambda*x:
        eigenvalue(mat A);
        // Define a generalized eigenvalue problem A*x = lambda*B*x:
        eigenvalue(mat A, mat B);
        
        // Define an eigenvalue problem of the form (M*lambda^2 + C*lambda + K)*u = 0:
        eigenvalue(mat K, mat C, mat M);
        // Define a general polynomial eigenvalue problem of the form
        // inmats[0] + inmats[1]*lambda + inmats[2]*lambda^2 + ...
        eigenvalue(std::vector<mat> inmats);
        
        void compute(int numeigenvaluestocompute, double targeteigenvaluemagnitude = 0.0);
        
        // Get the number of eigs found:
        int count(void);
        
        std::vector<double> geteigenvaluerealpart(void);
        std::vector<double> geteigenvalueimaginarypart(void);
        
        std::vector<vec> geteigenvectorrealpart(void);
        std::vector<vec> geteigenvectorimaginarypart(void);
        
        // Print the eigenvalues:
        void printeigenvalues(void);
        
        void printeigenfrequencies(void);
        
};

#endif
