// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <iostream>
#include <vector>
#include "polynomial.h"

class polynomial;

class polynomials
{
    private:

        std::vector<polynomial> mypolys = {};
        
        // Size in the ki, eta and phi direction:
        int mykilen = 0, myetalen = 0, myphilen = 0, mynummonomials = 0;
        std::vector<double> mycoeffs = {};

    public:
    
        polynomials(void) {};
        polynomials(std::vector<polynomial> input);
    
        int count(void) { return mypolys.size(); };
        
        // Evaluate at a single {ki,eta,phi} point:
        void evalatsingle(const std::vector<double>& evaluationpoint, std::vector<double>& evaled);
        
        
        void print(void);
        
};

#endif
