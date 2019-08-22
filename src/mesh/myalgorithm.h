// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef MYALGORITHM_H
#define MYALGORITHM_H

#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include "polynomial.h"

namespace myalgorithm
{
    // 'stablecoordinatesort' takes as input a vector of node coordinates 
    // in the format [coord1x coord1y coord1z coord2x ...] and provides a
    // vector that can be used to reorder the nodes so that the coordinates
    // are sorted according to first the x, then the y and finally the z 
    // coordinates: sortedcoordinates = coordinates(reorderingvector,:).
    // The sorting is not affected by roundoff noise up to a threshold.
    // The sorting algorithm is stable.
    void stablecoordinatesort(std::vector<double> noisethreshold, std::vector<double>& coordinates, std::vector<int>& reorderingvector);
    
    // 'removeduplicatedcoordinates' outputs a vector that can be used to 
    // renumber the nodes so that all duplicates are removed:
    // coordinateswithoutduplicates(renumberingvector,:) = coordinates.
    // The output gives the number of non-duplicated nodes.
    int removeduplicatedcoordinates(std::vector<double> noisethreshold, std::vector<double>& coordinates, std::vector<int>& renumberingvector);
    
    // 'stablesort' is like 'stablecoordinatesort' but instead of sorting 
    // according to coordinates it sorts according to 'tosort', a vector of ints.
    void stablesort(std::vector<int>& tosort, std::vector<int>& reorderingvector);
    
    // 'stablesortparallel' sorts in parallel according to tosort[0] then 
    // according to tosrt[1],... ALL int* arrays have size 'numentries'.
    // The output is a reordering vector.
    int* stablesortparallel(std::vector<int*> tosort, int numentries);
    
    std::vector<int> intersect(std::vector<int> a, std::vector<int> b);
    
    // Compressed sparse row to ijk format converter for sparse matrices.
    // 'csrrows' is in CSR format. The ijk row format is output in 'ijkrows'.
    // 'ijkrows' must be preallocated with a size equal to the nnz.
    void csrtoijk(int numberofrows, int* csrrows, int* ijkrows);
    
    // Get a single solution {ki,eta,phi} of a system of one up to three polynomials {poly1,poly2,...} equaled to the rhs.
    // The converged solution is placed in 'initialguess'. If at any iteration either ki, eta or phi
    // is farther away than 'boxsize' from the origin then the function stops and returns false (true otherwise).
    // The initial guess is supposed to be inside the box.
    bool getroot(std::vector<polynomial>& poly, std::vector<double>& rhs, std::vector<double>& initialguess, double boxsize = 3, double tol = 1e-12, int maxit = 50);
    
    // Split the 'tosplit' vector into 'blocklen' vectors of length tosplit.size()/blocklen.
    std::vector<std::vector<double>> splitvector(std::vector<double>& tosplit, int blocklen);
    
    // Norm each block in the vector:
    std::vector<double> normblocks(std::vector<double>& tonorm, int blocklen);
    
    // Get the file extension (dot included) in a string (works for any extension length):
    std::string getfileextension(std::string filename);
    // Get the file name without the path and the extension:
    std::string getfilename(std::string filename);
    
    // Return "s" for an argument bigger than 1, return "" otherwise:
    std::string getplurals(int count);
    
};

#endif
