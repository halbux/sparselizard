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

    public:

        polynomials(std::vector<polynomial> input);
        
};

#endif
