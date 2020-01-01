// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// Legendre polynomials are required to define the hierarchical 
// form functions used in Sabine Zaglmayr's thesis.

#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <iostream>
#include <vector>
#include "polynomial.h"

namespace legendre
{    
    // Legendre polynomials ln(x).
    // Output is ln for n = 0 ... maxn.
    std::vector<polynomial> l(int maxn, polynomial x);
    // Integrated Legendre polynomials Ln(x):
    std::vector<polynomial> L(int maxn, polynomial x);
    // Scaled Legendre polynomials ls(x,t):
    std::vector<polynomial> ls(int maxn, polynomial x, polynomial t);
    // Scaled integrated Legendre polynomials Ls(x,t):
    std::vector<polynomial> Ls(int maxn, polynomial x, polynomial t);
};

#endif
