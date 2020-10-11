// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This code calls BLAS. See http://www.netlib.org/blas/ for more information.

// This object stores the ROW-MAJOR values of a dense matrix of doubles.
// In other words the matrix is stored as [row1 row2 row3 ...].

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "intdensematrix.h"

class densematrix
{   

    private:

        long long int numrows = 0;
        long long int numcols = 0;

        // This is only used in 'multiply':
        bool istransposed = false;

        std::shared_ptr<double> myvalues = NULL;

        // Throws an error if matrix is empty:
        void errorifempty(void);

    public:

        // Set empty matrix:
        densematrix(void) {};
        // Set number of rows and columns:
        densematrix(long long int numberofrows, long long int numberofcolumns);
        // Initialise to a double value. Use this to set a matrix with constant value.
        densematrix(long long int numberofrows, long long int numberofcolumns, double initvalue);
        // Initialise with a vector (row major):
        densematrix(long long int numberofrows, long long int numberofcolumns, const std::vector<double> valvec);
        // Initialise to consecutive numbers [init init+step init+2*step ...].
        densematrix(long long int numberofrows, long long int numberofcolumns, double init, double step);
        // Vertical concatenation of dense matrices:
        densematrix(std::vector<densematrix> input);

        long long int countrows(void) { return numrows; };
        long long int countcolumns(void) { return numcols; };
        long long int count(void) { return numrows*numcols; };

        // Was it initialised?
        bool isdefined(void) { return (myvalues != NULL); };

        // Set the values of a given row.
        void setrow(long long int rownumber, std::vector<double> rowvals);

        // Output a flattened matrix, i.e. put all rows one after the other.
        // In practice this only requires to change 'numrows' and 'numcols'.
        // Values are NOT copied!
        densematrix flatten(void);

        // Insert a block in the matrix. Top left of block is at (row, col):
        void insert(long long int row, long long int col, densematrix toinsert);

        // Insert a matrix at given row numbers. The number of columns must be the same.
        void insertatrows(std::vector<int> selectedrows, densematrix toinsert);
        void insertatcolumns(std::vector<int> selectedcolumns, densematrix toinsert);

        // For fast access take into account that the storage is row-major.
        // Getting the values pointer with 'getvalues' is a better choice!
        double getvalue(long long int rownumber, long long int columnnumber);
        void setvalue(long long int rownumber, long long int columnnumber, double val);
        
        void getvalues(std::vector<double>& topopulate);

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
        // Add A*B without changing A or B.
        void addproduct(densematrix A, densematrix B);
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
        void tan(void);
        void asin(void);
        void acos(void);
        void atan(void);
        void log10(void);
        void mod(double modval);

        // Get the min and max values in out[0] and out[1] respectively:
        std::vector<double> minmax(void);

        // Get the max value:
        double max(void);
        // Get the max absolute value:
        double maxabs(void);

        // Sum all values:
        double sum(void);

        // Multiply column i by input[i].
        void multiplycolumns(std::vector<double> input);

        // [Arow1*Brow1  Arow1*Brow2  ... Arow2*Brow1...].
        densematrix multiplyallrows(densematrix input);

        // This special product is called by an el x (gp x ffd) matrix A where the columns are grouped by ffd blocks of gp columns.
        // The 'tfval' matrix has size fft x gp. The returned matrix has size el x (gp x ffd x fft) and corresponds to [A*tfvalrow1 A*tfvalrow2 ...].
        densematrix dofinterpoltimestf(densematrix tfval);
        
        // [A1 A2 ...].multiplycolumns(B) replaces the calling matrix by [A1*B A2*B ...] where Ai*B is the elementwise product of Ai and B.
        void multiplycolumns(densematrix input);
        
        // A becomes [A; A; A; ...] n times.
        densematrix duplicatevertically(int n);
        // A becomes [A A A ...] n times.
        densematrix duplicatehorizontally(int n);

        // Extract a set of rows/columns from the matrix:
        densematrix extractrows(std::vector<int> selected);
        densematrix extractcols(std::vector<int> selected);
        // Extract the rows/columns in a given range:
        densematrix extractrows(long long int rangebegin, long long int rangeend);
        densematrix extractcols(long long int rangebegin, long long int rangeend);
        
        // Multiply this block diagonal matrix (column major) by a vector:
        densematrix blockdiagonaltimesvector(intdensematrix blocklens, densematrix v);

};

#endif
