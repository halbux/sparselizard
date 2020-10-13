// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef MYALGORITHM_H
#define MYALGORITHM_H

#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "element.h"
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
    
    // 'removeduplicates' outputs a vector that can be used to 
    // renumber the nodes so that all duplicates are removed:
    // coordinateswithoutduplicates(renumberingvector,:) = coordinates.
    // The output gives the number of non-duplicated nodes.
    int removeduplicates(std::vector<double>& coordinates, std::vector<int>& renumberingvector);
    
    // This is for a vector of ints:
    void stablesort(std::vector<int>& tosort, std::vector<int>& reorderingvector);
    // This is for a vector of doubles:
    void stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector);
    // Same but sort by blocks of size 'blocklen':
    void stablesort(double noisethreshold, std::vector<double>& tosort, std::vector<int>& reorderingvector, int blocklen);
    
    void tuple3sort(std::vector<std::tuple<int,int,double>>& tosort);
    
    // Slice coordinates 'toslice' in the x, y and z direction into nsx*nsy*nsz groups. First slice position and distance between slices is provided as argument. Returned containers are:
    //
    //    - Group address 'ga[g]' gives the first position in 'pn' for group g (length of 'ga' is the number of groups + 1 and last value is the number of coordinates)  
    //    - Point number 'pn[i]' gives the corresponding point number in 'toslice' (must be preallocated to number of coordinates)
    //    - Point coordinates 'pc[3*i+0/+1/+2]' gives the x, y and z coordinates of point 'pn[i]' (must be preallocated to 3 x number of coordinates)
    //
    void slicecoordinates(std::vector<double>& toslice, double minx, double miny, double minz, double dx, double dy, double dz, int nsx, int nsy, int nsz, std::vector<int>& ga, int* pn, double* pc);
    
    // Return {xmin,xmax,ymin,ymax,zmin,zmax} for the coordinates {x1,y1,z1,x2,...}:
    std::vector<double> getcoordbounds(std::vector<double>& coordinates);
    
    std::vector<int> intersect(std::vector<int> a, std::vector<int> b);
    
    // Compressed sparse row to ijk format converter for sparse matrices.
    // 'csrrows' is in CSR format. The ijk row format is output in 'ijkrows'.
    // 'ijkrows' must be preallocated with a size equal to the nnz.
    void csrtoijk(int numberofrows, int* csrrows, int* ijkrows);
    
    // Get a single solution {ki,eta,phi} of a system of one up to three polynomials {poly1,poly2,...} equaled to the rhs.
    // After convergence 1 is returned and the solution is placed in 'initialguess'. If at any iteration either ki, eta or phi
    // is farther away than 'boxsize' from the origin then the function stops and returns 0. Value -1 is returned in any other case.
    // The initial guess is supposed to be inside the box.
    int getroot(polynomials& polys, std::vector<double>& rhs, std::vector<double>& initialguess, double boxsize = 2, double tol = 1e-10, int maxit = 20);

    // Get the reference coordinates and the element numbers corresponding to each (x,y,z) coordinate provided as argument. 
    // If the ith coordinate (xi,yi,zi) cannot be found in any element of the disjoint region then elems[i] is unchanged.
    // Any coordinate for which elems[i] is not -1 is ignored. 'elems' and 'kietaphis' must be preallocated to size numcoords and 3*numcoords.
    // This function is designed to be called in a for loop on multiple disjoint regions of same element type.
    void getreferencecoordinates(coordinategroup& coordgroup, int disjreg, std::vector<int>& elems, std::vector<double>& kietaphis);
 
    // Split the 'tosplit' vector into 'blocklen' vectors of length tosplit.size()/blocklen.
    std::vector<std::vector<double>> splitvector(std::vector<double>& tosplit, int blocklen);
    
    // Split a vector in two vectors according to 'select':
    void splitvector(std::vector<int>& vec, std::vector<bool>& select, std::vector<int>& falses, std::vector<int>& trues);
    
    // Norm each block in the vector:
    std::vector<double> normblocks(std::vector<double>& tonorm, int blocklen);
    
    // Return the interval number in which 'val' is. The interval tics must be sorted ascendingly.
    // Values out of bound are set to either the interval 0 or the last one. Size of 'tics' must be at least 2.
    int findinterval(double val, std::vector<double>& tics);
    // Get the tics for equally sized intervals:
    std::vector<double> getintervaltics(double minval, double maxval, int numintervals);
    
    // Get the file extension (dot included) in a string (works for any extension length):
    std::string getfileextension(std::string filename);
    // Get the file name without the path and the extension:
    std::string getfilename(std::string filename);
    
    // Return "s" for an argument bigger than 1, return "" otherwise:
    std::string getplurals(int count);
    
    // Get a vector of equally spaced numbers:
    std::vector<int> getequallyspaced(int start, int space, int amount);
    
    // Create a vector equal to n duplicates of the argument vector:
    std::vector<double> duplicate(std::vector<double> invec, int n);
    
    // Remove trailing CR (added in windows for newline):
    void osclean(std::string& line);
    
    // Convert an int to/from a boolean vector:
    std::vector<bool> inttobinary(int numbits, int num);
    int binarytoint(std::vector<bool> num);
    
    // Output a unique number for each combination of number inequalities (no equality allowed).
    // The output is a number ranging from 0 to numbers.size()! - 1.
    int identifyrelations(std::vector<int> numbers);
    
    // Factorial of positive integers or zero:
    int factorial(int n);
    
    // Provide the corner coordinates of each element concatenated in each element type.
    // This function returns (flattened from lowest type to highest) the edge number of
    // each edge in an element as well as a bool whose value is true if the edge barycenter
    // is close enough to any node in the corner coordinates.
    void assignedgenumbers(std::vector<std::vector<double>>& cornercoords, std::vector<int>& edgenumbers, std::vector<bool>& isbarycenteronnode);
    
    // For a vector 'vec' of repeating blocks [b0 b1 b2 ...] the output is [b0[sel[0]] b1[sel[0]] ... b0[sel[1]] b1[sel[1]] ...].
    std::vector<double> separate(std::vector<double>& v, int blocklen, std::vector<int> sel);
    
    // Return the renumbering obtained when the original renumbering is renumbered again by a new renumbering:
    std::vector<int> chainrenumbering(std::vector<int>& originalrenum, std::vector<int>& newrenum);
    
    // Invert the renumbering:
    std::vector<int> invertrenumbering(std::vector<int>& renum);
    // Get the reordering vector corresponding to the renumbering vector:
    std::vector<int> getreordering(std::vector<int>& renum);
    
    // Reorder an address-data pair based on a renumbering of the addresses (last entry in address vectors is data.size):
    void reorder(std::vector<int>& inad, std::vector<double>& indat, std::vector<int>& renumbering, std::vector<int>& outad, std::vector<double>& outdat);
    
    // 'elems' has format {elemtype0,elemnum0,elemtype1,...}. 'totalnumelems[i]' gives the number of elements of type i in the mesh.
    // For the ith point in 'elems' the vector indexinrcsoforigin[i]' gives the index/3 in 'rcs[element type at ith point in elems]'. 
    void toaddressdata(std::vector<int>& elems, std::vector<double>& refcoords, std::vector<int> totalnumelems, std::vector<std::vector<int>>& ads, std::vector<std::vector<double>>& rcs, std::vector<int>& indexinrcsoforigin);
    
    // Concatenate vectors:
    std::vector<int> concatenate(std::vector<std::vector<int>> tocat);
    
    // Return -1 if a < b, 0 if a = b and +1 if a > b:
    int inequalitytoint(int a, int b);
    
    // Norm a vector:
    void normvector(std::vector<double>& tonorm);
    
};

#endif
