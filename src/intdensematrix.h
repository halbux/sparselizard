// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object stores the ROW-MAJOR values of a dense matrix of int.
// In other words the matrix is stored as [row1 row2 row3 ...].

#ifndef INTDENSEMATRIX_H
#define INTDENSEMATRIX_H

#include <iostream>
#include <vector>
#include <memory>

class intdensematrix
{   

    private:
        
        long long int numrows = 0;
        long long int numcols = 0;
        
        std::shared_ptr<int> myvalues = NULL;
        
    public:
        
        // Set empty matrix:
        intdensematrix(void) {};
        // Set number of rows and columns:
        intdensematrix(long long int numberofrows, long long int numberofcolumns);
        // Initialise to a value:
        intdensematrix(long long int numberofrows, long long int numberofcolumns, int initvalue);
        // Initialise with a vector (row major):
        intdensematrix(long long int numberofrows, long long int numberofcolumns, const std::vector<int> valvec);
        // Initialise to consecutive numbers [init init+step init+2*step ...].
        intdensematrix(long long int numberofrows, long long int numberofcolumns, int init, int step);
        
        long long int countrows(void) { return numrows; };
        long long int countcolumns(void) { return numcols; };
        long long int count(void) { return numrows*numcols; };
        
        // Count the number of positive or zero integer values:
        long long int countpositive(void);

        void print(void);
        void printsize(void);
 
        int* getvalues(void);

        // A.duplicateallrowstogether(int n) for a matrix A of size pxq 
        // whose form is  [row1; row2; row3; ...] outputs a matrix of size 
        // (p*n)xq where every row of matrix A has been duplicated n 
        // times as follows [row1; row2; row3;  ... row1; row2; row3; ...].
        intdensematrix duplicateallrowstogether(int n);
        // A.duplicaterowsonebyone(int n) for a matrix A of size pxq 
        // whose form is [row1; row2; row3; ...] outputs a matrix of size 
        // (p*n)xq where every row of matrix A has been duplicated n 
        // times as follows [row1; row1; row1;  ... row2; row2; row2; ...].
        intdensematrix duplicaterowsonebyone(int n);

};

#endif
