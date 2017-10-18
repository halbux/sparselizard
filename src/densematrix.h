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
        
        // Cannot be used with all functions:
        bool istransposed = false;
        
        std::shared_ptr<double> myvalues = NULL;
        
	public:
        
        // Set empty matrix:
        densematrix(void) {};
        // Set number of rows and columns:
        densematrix(int numberofrows, int numberofcolumns);
        // Initialise to a double value. Use this to set a matrix with constant value.
        densematrix(int numberofrows, int numberofcolumns, double initvalue);
        // Initialise to a row or column vector:
        densematrix(int numberofrows, int numberofcolumns, std::vector<double>& valvec);
        // Initialise to consecutive numbers [init init+step init+2*step ...].
        // Usefull for debug purposes at least.
        densematrix(int numberofrows, int numberofcolumns, int init, int step);
        
        int countrows(void) { return numrows; };
        int countcolumns(void) { return numcols; };

        // Was it initialised?
        bool isdefined(void) { return (numrows != 0 && numcols != 0); };
        
        // Set a given row number to a value.
        void setrow(int, std::vector<double>); // REMOVE!

        // Insert a block in the matrix. Top left of block is at (row, col):
        void insert(int row, int col, densematrix toinsert);
        
        // Insert a matrix at given row numbers. The number of columns must be the same.
        void insertatrows(std::vector<int> selectedrows, densematrix toinsert);
        void insertatcolumns(std::vector<int> selectedcolumns, densematrix toinsert);
        
        // Output a flattened matrix, i.e. put all rows one after the other.
        // In practice this only requires to change 'numrows' and 'numcols'.
        densematrix flatten(void);

        // For fast access take into account that the storage is 
        // row-major for the not-transposed case.
        double getvalue(int rownumber, int columnnumber);
        void setvalue(int rownumber, int columnnumber, double);
        
//         bool operator==(densematrix input);
        
        // Get a full copy, including the pointed values.
        densematrix copy(void);
        void print(void);
        void printsize(void);
        void transpose(void);
        
        bool isallzero(void);
        
        // Multiply current object matrix by B:
        densematrix multiply(densematrix);
        
        // The matrix cannot get out of scope
        double* getvalues(void);
        
        densematrix getcolumns(std::vector<int> cols);
        
        // Elementwise operations. They do not accept matrix transposition.
        // Matrices must all have the same size.
        void erroriftransposed(void);
        void multiplyelementwise(densematrix);
        void add(densematrix);
        // Multiply by coef*B without changing B.
        void addproduct(double coef, densematrix B);
        densematrix returnproduct(double coef);
        void subtractelementwise(densematrix);
        void multiplyelementwise(double);
        void minus(void);
        void power(densematrix);
        // Compute 1/val for all values:
        void invert(void);
        void abs(void);
        void sin(void);
        void cos(void);
        void log10(void);
        
        double maxabs(void);

        double sum(void);
        
        void multiplycolumns(std::vector<double>);
        
        // [Acol1*Bcol1  Acol1*Bcol2  ... Acol2*Bcol1...]
        densematrix multiplyallrows(densematrix);
        
        // mat becomes [mat; mat; mat; ...] n times
        densematrix duplicatevertically(int n);

};

#endif
