// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

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
#include <slepceps.h>

class eigenvalue
{
    private:
        
        mat myA;
        mat myB;
        
        // Real and imaginary part of the eigenvalues and eigenvectors:
        std::vector<double> eigenvaluereal = {};
        std::vector<double> eigenvalueimaginary = {};
        
        std::vector<vec> eigenvectorreal = {};
        std::vector<vec> eigenvectorimaginary = {};
        
    public:

        // Define a standard eigenvalue problem A*x = lambda*x:
        eigenvalue(mat A);
        // Define a generalised eigenvalue problem A*x = lambda*B*x:
        eigenvalue(mat A, mat B);
        
        void compute(int numeigenvaluestocompute, double targeteigenvaluemagnitude = 0.0);
        
        // Get the number of eigs found:
        int count(void) { return eigenvaluereal.size(); };
        
        std::vector<double> geteigenvaluerealpart();
        std::vector<double> geteigenvalueimaginarypart();
        
        std::vector<vec> geteigenvectorrealpart();
        std::vector<vec> geteigenvectorimaginarypart();
        
        // Print the eigenvalues:
        void printeigenvalues(void);
        // In case the eigenvalues are real and correspond to the square of the angular frequency:
        void printeigenfrequencies(void);
        
};

#endif
