// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef GENTOOLS_H
#define GENTOOLS_H

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
#include "densemat.h"

namespace gentools
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
    
    // Provide a renumbering vector to remove duplicated coordinates.
    // Entry renumberingvector[i] is guaranteed lower or equal to i.
    // The output gives the number of non-duplicated coordinates.
    int removeduplicates(std::vector<double>& coordinates, std::vector<int>& renumberingvector);
    // Provide a renumbering vector to remove duplicated blocks.
    // Entry renumberingvector[i] is guaranteed lower or equal to i.
    // The output gives the number of non-duplicated blocks.
    int removeduplicates(std::vector<int> toremove, std::vector<int>& renumberingvector, int blocklen);
    
    // Remove duplicated coordinates:
    void removeduplicates(std::vector<double>& coordinates);
    
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
    
    std::vector<int> unique(std::vector<int> a);
    
    std::vector<int> intersect(std::vector<int> a, std::vector<int> b);
    std::vector<int> exclude(std::vector<int> a, std::vector<int> b);
    
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
    
    // Select indexes of a vector:
    void select(std::vector<int>& vals, std::vector<int>& selectedindexes, std::vector<int>& selected);
    void select(std::vector<bool>& vals, indexmat selectedindexes, std::vector<bool>& selected);
    
    // Compare the ordering of two vectors (vectors must be circularly or anti-circularly identical and not empty).
    // Length 1 is considered not flipped. Length 2 is considered flipped if not identical.
    bool isflipped(std::vector<int>& a, std::vector<int>& b);
    
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
    
    // Return the exact conversion from integer to double (throw an error if not possible):
    double exactinttodouble(long long int num);
    
    // Output a unique number for each combination of number inequalities (no equality allowed).
    // The output is a number ranging from 0 to numbers.size()! - 1.
    int identifyrelations(std::vector<int> numbers);
    
    // Factorial of positive integers or zero:
    int factorial(int n);
    
    // Provide the corner coordinates of each element concatenated in each element type.
    // This function returns (flattened from lowest type to highest) the edge number of
    // each edge in an element as well as a bool whose value is true if the edge barycenter
    // is close enough to any node in the corner coordinates. In the DDM framework the own
    // elements must include all / the inner overlap elements for no-overlap / overlap DDM.
    void assignedgenumbers(std::vector<bool>& isownelem, std::vector<std::vector<double>>& cornercoords, std::vector<int>& edgenumbers, std::vector<bool>& isbarycenteronnode);
    
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
    void concatenate(std::vector<std::vector<double>>& tocat, std::vector<double>& cated);
    
    // Return -1 if a < b, 0 if a = b and +1 if a > b:
    int inequalitytoint(int a, int b);
    
    // Norm a vector:
    void normvector(std::vector<double>& tonorm);
    
    // Solve Ux = b where U is column-major upper triangular {r0c0,r0c1,r1c1,r0c2,...}:
    void solveuppertriangular(int len, double* U, double* b, double* x);

    // Givens rotation:
    void givensrotation(double a, double b, double& c, double& s, double& r);
    
    // Apply Givens rotation to h (the kth column of the unreduced upper Hessenberg matrix, h must have length k+2):
    void applygivensrotation(double* h, std::vector<double>& cs, std::vector<double>& sn, int k);
    
    // Q has one row per Krylov vector (at least k+2 rows must be preallocated).
    // Column k of the unreduced upper Hessenberg matrix is returned (length k+2):
    std::vector<double> arnoldi(densemat (*mymatmult)(densemat), densemat Q, int k);
    
    // Create a vector to renumber integer values with gaps to values without gap.
    // Integers must all be positive or zero. The number of unique integers is returned.
    int squeeze(std::vector<int>& nums, int maxval, std::vector<int>& renumbering);
    
    // Set 'numtrue' to -1 if the number of true entries is not known:
    void find(std::vector<bool>& invec, int numtrue, std::vector<int>& trueindexes);
    
    // For each coordinate to find this gives the index of the coordinate matched in the target (-1 if not found in the target).
    // The number of matches is returned. It is allowed to have duplicate coordinates in the target.
    int findcoordinates(std::vector<double>& targetcoords, std::vector<double>& tofindintarget, std::vector<int>& posfound);
    
    // Write to 'selectedcoords' the coordinates in 'coords' that have a true selection value: 
    void selectcoordinates(std::vector<bool>& selection, std::vector<double>& coords, double* selectedcoords);
    
    // From the candidates vector pick a set of coordinates that have rather equally spaced (possibly NOT UNIQUE) indexes:
    void pickcandidates(int numbertopick, std::vector<double>& candidatecoordinates, std::vector<double>& picked);
    
    // 'data[i][j]' is a vector that can either have size 0 or a size larger or equal to 2. In the latter case the vector
    // is the result of the concatenation of two vectors a and b and has format {lengtha, lengthb, ... valuesa ..., ... valuesb ...}.
    // This function extracts the sizes of each a/b vector into 'sizesa'/'sizesb' and concatenates 
    // all valuesa/valuesb vectors together into the output 'dataa'/'datab'.
    void split(std::vector< std::vector<std::vector<double>> >& data, std::vector<double>& dataa, std::vector<double>& datab, std::vector<std::vector<int>>& sizesa, std::vector<std::vector<int>>& sizesb);
    
    // This function creates 'packed' of size numtags where 'packed[k]' holds all data associated to the kth tag.
    // 'topack' has size numtags^2 and 'topack[i*numtags+j]' (only j >= i+1 is considered, empty vectors are skipped)
    // holds data of tag 'tags[j]' associated to the ith tag and data of tag 'tags[i]' associated to the jth tag.
    // 'packed[k]' is the concatenation of all data from 'topack' associated to the kth tag. It has format
    // {tag0, length0, ... data0 ..., tag1, length1, ... data1 ..., ...}.
    void pack(std::vector<int> tags, std::vector<std::vector<double>>& topack, std::vector<std::vector<double>>& packed);
    
    // Unpack 'packed' of format {tag0, length0, ... data0 ..., tag1, length1, ... data1 ..., ...} into 'unpack'.
    // 'unpacked' has a size equal to the number of tags (possibly not unique) found in 'packed'.
    // 'unpacked[i]' holds the ith dataset (of tag 'output[i]') found in 'packed'.
    // All tags found are returned in the order of appearance. They might not be unique.  
    std::vector<int> unpack(std::vector<double>& packed, std::vector<std::vector<double>>& unpacked);
    
    // Return all values at data[i*period+shift] and remove them from 'data' (length of 'data' should be a multiple of the period):
    std::vector<int> extract(std::vector<int>& data, int period, int shift);
    std::vector<double> extract(std::vector<double>& data, int period, int shift);
    
    // Return the ceiled division a/b:
    int ceildiv(int a, int b);
    
    // Get the number of integers in the packed vector:
    int getpackedsize(int numbits);
    // Pack a vector of booleans (works for an arbitrary sizeof(int) bytes):
    void pack(std::vector<bool>& topack, std::vector<int>& packed);
    // Unpack the vector of booleans (was of length 'orignumbools' before being packed):
    void unpack(int orignumbools, std::vector<int>& packed, std::vector<bool>& unpacked);
    
    // Split values < val and >= val while keeping the ordering:
    void split(std::vector<int>& vals, int val, std::vector<int>& lowervals, std::vector<int>& highervals);
    
    // Count the number of true entries:
    int counttrue(std::vector<bool>& tocount);
    
    // For a vector of positive or zero values compress consecutive zero values:
    void compresszeros(std::vector<int>& tocompress);
    void decompresszeros(std::vector<int>& todecompress);
    
    // Get the maximum element dimension found in the element list (-1 if empty):
    int getmaxdim(std::vector<std::vector<int>>* elementlist);
    
    // Sum all entries in the vector (return 0 if empty):
    int sum(std::vector<int>& tosum);
    double sum(std::vector<double>& tosum);
    
    // Split a string at the first colon. If there is no colon then the string to split is placed in 'last':
    void splitatcolon(std::string tosplit, std::string& first, std::string& last);
    
    // Find the true and false indexes in the argument vector and provide the renumbering of each vector entry to its true/false index:
    void findtruefalse(std::vector<bool>& invec, indexmat& trueinds, indexmat& falseinds, std::vector<int>& renum);
    
    // For every edge in 'physreg' (in the 'der' order) return a flag telling if its direction has to be flipped
    // to fullfill the condition that at every node the touching edges point together either inwards or outwards.
    void inoutorient(int physreg, std::vector<bool>& flipit);
    // Helper function to be called recursively.
    void inoutorient(int startnode, std::vector<int>& edgestatus, bool isoutward, bool isrecursivecall);
    
    // The inner overlap cell values are decided by the owner of the inner overlap:
    void fixatoverlap(std::vector<std::vector<int>>& cellvalues);
    
    // Get the list of edges in the inner overlap interface with each neighbour and preallocate for the outer overlap interface.
    // In case of no-overlap DDM the inner and outer overlap interfaces both equal the no-overlap interface.
    void getedgesininnerinterfaces(std::vector<std::vector<int>>& iiedgelists, std::vector<std::vector<int>>& oiedgelistspreallocated);
    
    // Append the values of this rank and all neighbour ranks (in case of DDM). Also get the value 'togroup' on each neighbour rank.
    std::vector<int> appendneighbourvalues(std::vector<double>& toappendto, std::vector<double>& toappend, int togroup);
    
    // From the provided disjoint regions (with same element type) return the elements with at least one active corner node:
    std::vector<int> getactiveelements(std::vector<int> disjregs, std::vector<bool>& isnodeactive);
};

#endif
