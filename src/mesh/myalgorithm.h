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
#include "polynomials.h"
#include "coordinategroup.h"

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
    // Same as above but first sort according to a vector of integers:
    void stablecoordinatesort(std::vector<double> noisethreshold, std::vector<int>& elems, std::vector<double>& coordinates, std::vector<int>& reorderingvector);
    
    // 'removeduplicatedcoordinates' outputs a vector that can be used to 
    // renumber the nodes so that all duplicates are removed:
    // coordinateswithoutduplicates(renumberingvector,:) = coordinates.
    // The output gives the number of non-duplicated nodes.
    int removeduplicatedcoordinates(std::vector<double> noisethreshold, std::vector<double>& coordinates, std::vector<int>& renumberingvector);
    
    // This is for a vector of ints:
    void stablesort(std::vector<int>& tosort, std::vector<int>& reorderingvector);
    // This is for a vector of doubles:
    void stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector);
    // Same but sort by blocks of size 'blocklen':
    void stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector, int blocklen);
    
    // 'stablesortparallel' sorts in parallel according to tosort[0] then 
    // according to tosrt[1],... ALL int* arrays have size 'numentries'.
    // The output is a reordering vector.
    int* stablesortparallel(std::vector<int*> tosort, int numentries);
    void tuple3sort(std::vector<std::tuple<int,int,double>>& tosort);
    
    // The output 'slices[i]' gives the indexes of all values of 'toslice' that are in the interval between minval+i*delta and minval+(i+1)*delta.
    // Value minval must be smaller than the smallest value in 'toslice' and minval+numslices*delta must be larger than the largest value.
    void slicecoordinates(double noisethreshold, std::vector<double>& toslice, double minval, double delta, int numslices, std::vector<std::vector<int>>& slices);
    
    // Return {xmin,xmax,ymin,ymax,zmin,zmax} for the coordinates {x1,y1,z1,x2,...}:
    std::vector<double> getcoordbounds(std::vector<double>& coordinates);
    
    std::vector<int> intersect(std::vector<int> a, std::vector<int> b);
    
    // Compressed sparse row to ijk format converter for sparse matrices.
    // 'csrrows' is in CSR format. The ijk row format is output in 'ijkrows'.
    // 'ijkrows' must be preallocated with a size equal to the nnz.
    void csrtoijk(int numberofrows, int* csrrows, int* ijkrows);
    
    // Get a single solution {ki,eta,phi} of a system of one up to three polynomials {poly1,poly2,...} equaled to the rhs.
    // The polynomials in 'polys' must be ordered as {poly1,dkipoly1,detapoly1,...,poly2,dkipoly2,...}.
    // After convergence 1 is returned and the solution is placed in 'initialguess'. If at any iteration either ki, eta or phi
    // is farther away than 'boxsize' from the origin then the function stops and returns 0. Value -1 is returned in any other case.
    // The initial guess is supposed to be inside the box.
    int getroot(std::vector<polynomial>& poly, std::vector<double>& rhs, std::vector<double>& initialguess, double boxsize = 2, double tol = 1e-12, int maxit = 20);
    int getroot(polynomials& polys, std::vector<double>& rhs, std::vector<double>& initialguess, double boxsize = 2, double tol = 1e-12, int maxit = 20);

    // Attempt to get the above mentionned root with multiple initial guesses provided in format {ki1,eta1,phi1,ki2,eta2,phi2,...}.
    int getrootmultiguess(std::vector<polynomial>& poly, std::vector<double>& rhs, std::vector<double>& initialguesses, std::vector<double>& kietaphi, double boxsize = 2, double tol = 1e-12, int maxit = 20);

    // Get the reference coordinates and the element numbers corresponding to each (x,y,z) coordinate provided as argument. 
    // If the ith coordinate (xi,yi,zi) cannot be found in any element of the disjoint region then elems[i] is unchanged.
    // Any coordinate for which elems[i] is not -1 is ignored. 'elems' and 'kietaphis' must be preallocated to size numcoords and 3*numcoords.
    // This function is designed to be called in a for loop on multiple disjoint regions of same element type.
    void getreferencecoordinates(coordinategroup& coordgroup, int disjreg, std::vector<int>& elems, std::vector<double>& kietaphis);
 
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
