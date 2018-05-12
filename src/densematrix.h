// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.

// This code calls BLAS. See http://www.netlib.org/blas/ for more information.

// This object stores the ROW-MAJOR values of a dense matrix of doubles.
// In other words the matrix is stored as [row1 row2 row3 ...].

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include "cblas.h"
#include <memory>

class densematrix
{   

    private:
        
        int numrows = 0;
        int numcols = 0;
        
        // This is only used in 'multiply':
        bool istransposed = false;
        
        std::shared_ptr<double> myvalues = NULL;

        // Throws an error if matrix is empty:
        void errorifempty(void);
        
    public:
        
        // Set empty matrix:
        densematrix(void) {};
        // Set number of rows and columns:
        densematrix(int numberofrows, int numberofcolumns);
        // Initialise to a double value. Use this to set a matrix with constant value.
        densematrix(int numberofrows, int numberofcolumns, double initvalue);
        // Initialise with a vector (row major):
        densematrix(int numberofrows, int numberofcolumns, const std::vector<double> valvec);
        // Initialise to consecutive numbers [init init+step init+2*step ...].
        densematrix(int numberofrows, int numberofcolumns, int init, int step);
        
        int countrows(void) { return numrows; };
        int countcolumns(void) { return numcols; };

        // Was it initialised?
        bool isdefined(void) { return (myvalues != NULL); };
        
        // Set the values of a given row.
        void setrow(int rownumber, std::vector<double>);
        
        // Output a flattened matrix, i.e. put all rows one after the other.
        // In practice this only requires to change 'numrows' and 'numcols'.
        // Values are NOT copied!
        densematrix flatten(void);

        // Insert a block in the matrix. Top left of block is at (row, col):
        void insert(int row, int col, densematrix toinsert);
        
        // Insert a matrix at given row numbers. The number of columns must be the same.
        void insertatrows(std::vector<int> selectedrows, densematrix toinsert);
        void insertatcolumns(std::vector<int> selectedcolumns, densematrix toinsert);

        // For fast access take into account that the storage is row-major.
        // Getting the values pointer with 'getvalues' is a better choice!
        double getvalue(int rownumber, int columnnumber);
        void setvalue(int rownumber, int columnnumber, double val);
        
        // Get a full copy (all values are copied).
        densematrix copy(void);
        
        void print(void);
        void printsize(void);
        
        bool isallzero(void);

        // Transpose has an effect only on the 'multiply' function below.
        // It does not actually transpose the matrix.
        void transpose(void);

        // Multiply current object matrix by B with BLAS:
        densematrix multiply(densematrix B);
        
        // The matrix cannot get out of scope
        double* getvalues(void);
        
        // Add coef*B without changing B.
        void addproduct(double coef, densematrix B);
        // Return the product by 'coef' without modifying this object.
        densematrix returnproduct(double coef);
        
        // Elementwise operations below. 
        // They do not accept matrix transposition.
        // Matrices must all have the same size.
        void multiplyelementwise(densematrix B);
        void multiplyelementwise(double val);
        void add(densematrix B);
        void subtract(densematrix B);
        void minus(void);
        void power(densematrix exponent);
        // Compute 1/val for all values:
        void invert(void);
        void abs(void);
        void sin(void);
        void cos(void);
        void log10(void);
        
        // Get the max value:
        double max(void);
        // Get the index of the max value:
        int maxindex(void);

        // Get the max absolute value:
        double maxabs(void);

        // Sum all values:
        double sum(void);
        
        // Multiply column i by input[i].
        void multiplycolumns(std::vector<double> input);
        
        // [Arow1*Brow1  Arow1*Brow2  ... Arow2*Brow1...].
        densematrix multiplyallrows(densematrix input);
        
        // A becomes [A; A; A; ...] n times.
        densematrix duplicatevertically(int n);

};

#endif
