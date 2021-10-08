// sparselizard - Copyright (C) see copyright file.
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

        int mynumpolys = 0;
        
        // Size in the ki, eta and phi direction:
        int mykilen = 0, myetalen = 0, myphilen = 0, mynummonomials = 0;
        std::vector<double> mycoeffs = {};

    public:
    
        polynomials(void) {};
        polynomials(std::vector<polynomial> input);
    
        int count(void) { return mynumpolys; };
        
        // Evaluate at a single {ki,eta,phi} point:
        void evalatsingle(const std::vector<double>& evaluationpoint, std::vector<double>& evaled);
        // Same as above but {poly1,dkipoly1,detapoly1,...,poly2,dkipoly2,detapoly2,...} is returned. 
        // 'num' equal to 0/1/2/3 returns respectively poly/poly+dki/poly+dki+deta/poly+dki+deta+dphi. 
        void evalatsingle(const std::vector<double>& evaluationpoint, int num, std::vector<double>& evaled);
        
        // Return the weighted sum of the original polynomials.
        // Output holds {p1,p2,...} where pi = sum_k( weights[ i * mynumpolys + k ] * originalpoly[k] ).
        polynomials sum(std::vector<double>& weights);
        
        void print(void);
        
};

#endif
