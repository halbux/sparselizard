// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <vector>
#include "polynomials.h"

class polynomial
{
    friend class polynomials;
    
    private:

        // 'mycoefficients[i][j][k]' gives the coefficient of the monomial ki^i*eta^j*phi^k
        std::vector<std::vector<std::vector<double>>> mycoefficients = {};

    public:

        void set(const std::vector<std::vector<std::vector<double>>>& coefficients);

        std::vector<std::vector<std::vector<double>>> get(void) { return mycoefficients; };

        // 'print' prints the polynomial to the console:
        void print(void);
        // Compute for the ki, eta and phi values provided.
        // Format is [ki1 eta1 phi1 ki2 eta2 phi2 ...].
        // Set the int to 0 to get the no derivative value, 1 for dki, 2 for deta and 3 for dphi.
        std::vector<double> evalat(const std::vector<double>& evaluationpoints, int whichderivative);
        // Defining the +, - and * operators for polynomials:
        polynomial operator*(polynomial);
        polynomial operator+(polynomial);
        polynomial operator-(polynomial);
        polynomial operator+();
        polynomial operator-();
        // Define the right product poly*double:
        polynomial operator*(double);
        polynomial operator+(double);
        polynomial operator-(double);

        // Defining the ki, eta and phi derivatives:
        void dki(void);
        void deta(void);
        void dphi(void);
        // Defining a generic derivative operator.
        // Argument set to 0 gives dki, 1 deta and 2 dphi.
        polynomial derivative(int whichderivative);
};

// Define the left version of the operators based on the right one.
polynomial operator*(double, polynomial);
polynomial operator+(double, polynomial);
polynomial operator-(double, polynomial);

#endif
