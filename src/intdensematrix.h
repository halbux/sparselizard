// sparselizard - Copyright (C) see copyright file.
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
        
        // Throws an error if matrix is empty:
        void errorifempty(void);
        
    public:
        
        // Set empty matrix:
        intdensematrix(void) {};
        // Set number of rows and columns:
        intdensematrix(long long int numberofrows, long long int numberofcolumns);
        // Initialise to a value:
        intdensematrix(long long int numberofrows, long long int numberofcolumns, int initvalue);
        // Initialise with a vector (row major):
        intdensematrix(long long int numberofrows, long long int numberofcolumns, std::vector<int> valvec);
        // Initialise to consecutive numbers [init init+step init+2*step ...].
        intdensematrix(long long int numberofrows, long long int numberofcolumns, int init, int step);
        // Vertical concatenation of dense matrices:
        intdensematrix(std::vector<intdensematrix> input);
        
        long long int countrows(void) { return numrows; };
        long long int countcolumns(void) { return numcols; };
        long long int count(void) { return numrows*numcols; };
        
        // Output the mxn resized matrix (this only changes 'numrows' and 'numcols'). Values are NOT copied!
        intdensematrix getresized(long long int m, long long int n);
        
        // Count the number of positive or zero integer values:
        long long int countpositive(void);
        // Count the number of occurences of a value:
        long long int countoccurences(long long int value);
        // Filter out the argument value and return a column vector:
        intdensematrix removevalue(long long int toremove);
        
        // Return a vector whose ith entry gives the number of times value i appears:
        std::vector<int> countalloccurences(int maxintval);
        // And all indexes at which it appears:
        std::vector<std::vector<int>> findalloccurences(int maxintval);

        // Sum all values:
        long long int sum(void);
        
        // Get the min and max values in out[0] and out[1] respectively:
        std::vector<int> minmax(void);
        
        int max(void);

        void print(void);
        void printsize(void);
 
        int* getvalues(void);
        
        // Get a full copy (all values are copied).
        intdensematrix copy(void);
        
        // Get the transpose without modifying this object.
        intdensematrix gettranspose(void);

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
        // Same for columns:
        intdensematrix duplicateallcolstogether(int n);
        intdensematrix duplicatecolsonebyone(int n);

        // Extract a set of rows/columns from the matrix:
        intdensematrix extractrows(std::vector<int>& selected);
        intdensematrix extractcols(std::vector<int>& selected);
        
        // Select all indexes for which the selection equals the last argument:
        intdensematrix select(std::vector<bool>& sel, bool selectif);

};

#endif
