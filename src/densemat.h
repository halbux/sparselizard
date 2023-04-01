// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This code calls BLAS. See http://www.netlib.org/blas/ for more information.

// This object stores the ROW-MAJOR values of a dense matrix of doubles.
// In other words the matrix is stored as [row1 row2 row3 ...].

#ifndef DENSEMAT_H
#define DENSEMAT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "logs.h"
#include "petscmat.h"
#include "indexmat.h"

class densemat
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
        densemat(void) {};
        // Set number of rows and columns:
        densemat(long long int numberofrows, long long int numberofcolumns);
        // Initialise to a double value. Use this to set a matrix with constant value.
        densemat(long long int numberofrows, long long int numberofcolumns, double initvalue);
        // Initialise with a vector (row major):
        densemat(long long int numberofrows, long long int numberofcolumns, std::vector<double> valvec);
        // Initialise to consecutive numbers [init init+step init+2*step ...].
        densemat(long long int numberofrows, long long int numberofcolumns, double init, double step);
        // Vertical concatenation of dense matrices:
        densemat(std::vector<densemat> input);

        long long int countrows(void) { return numrows; };
        long long int countcolumns(void) { return numcols; };
        long long int count(void) { return numrows*numcols; };

        // Was it initialised?
        bool isdefined(void) { return (myvalues != NULL); };

        // Set the values of a given row.
        void setrow(long long int rownumber, std::vector<double> rowvals);
        
        // Output the mxn resized matrix (this only changes 'numrows' and 'numcols'). Values are NOT copied!
        densemat getresized(long long int m, long long int n);
        // Output the 1x(m*n) resized matrix:
        densemat getflattened(void);

        // Insert a block in the matrix. Top left of block is at (row, col):
        void insert(long long int row, long long int col, densemat toinsert);

        // Insert a matrix at given row numbers. The number of columns must be the same.
        void insertatrows(std::vector<int> selectedrows, densemat toinsert);
        void insertatcolumns(std::vector<int> selectedcolumns, densemat toinsert);

        // For fast access take into account that the storage is row-major.
        // Getting the values pointer with 'getvalues' is a better choice!
        double getvalue(long long int rownumber, long long int columnnumber);
        void setvalue(long long int rownumber, long long int columnnumber, double val);
        
        void getvalues(std::vector<double>& topopulate);

        // Get a full copy (all values are copied).
        densemat copy(void);

        void print(void);
        void printsize(void);

        bool isallzero(void);

        // Transpose has an effect only on the 'multiply' function below.
        // It does not actually transpose the matrix.
        void transpose(void);

        // Multiply current object matrix by B with BLAS:
        densemat multiply(densemat B);

        // The matrix cannot get out of scope
        double* getvalues(void);

        // Add coef*B without changing B.
        void addproduct(double coef, densemat B);
        // Add A*B without changing A or B.
        void addproduct(densemat A, densemat B);
        // Get the product by 'coef' without modifying this object.
        densemat getproduct(double coef);

        // Get the transpose without modifying this object.
        densemat gettranspose(void);
        
        // Get the matrix inverse (must be square):
        densemat getinverse(void);

        // Elementwise operations below. 
        // Matrices must all have the same size.
        void multiplyelementwise(densemat B);
        void multiplyelementwise(double val);
        void add(densemat B);
        void subtract(densemat B);
        void minus(void);
        void power(densemat exponent);
        // Compute 1/val for all values:
        void invert(void);
        void abs(void);
        void sin(void);
        void cos(void);
        void tan(void);
        void asin(void);
        void acos(void);
        void atan(void);
        void log(void);
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
        densemat multiplyallrows(densemat input);

        // This special product is called by an el x (gp x ffd) matrix A where the columns are grouped by ffd blocks of gp columns.
        // The 'tfval' matrix has size fft x gp. The returned matrix has size el x (gp x ffd x fft) and corresponds to [A*tfvalrow1 A*tfvalrow2 ...].
        densemat dofinterpoltimestf(densemat tfval);
        
        // [A1 A2 ...].multiplycolumns(B) replaces the calling matrix by [A1*B A2*B ...] where Ai*B is the elementwise product of Ai and B.
        void multiplycolumns(densemat input);
        
        // A becomes [A; A; A; ...] n times.
        densemat duplicatevertically(int n);
        // A becomes [A A A ...] n times.
        densemat duplicatehorizontally(int n);

        // Extract a set of rows/columns from the matrix:
        densemat extractrows(std::vector<int>& selected);
        densemat extractcols(std::vector<int>& selected);
        // Extract the rows/columns in a given range:
        densemat extractrows(long long int rangebegin, long long int rangeend);
        densemat extractcols(long long int rangebegin, long long int rangeend);
        
        // Multiply this block diagonal matrix (column major) by a vector:
        densemat blockdiagonaltimesvector(indexmat blocklens, densemat v);

};

#endif
