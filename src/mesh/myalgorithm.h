// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef MYALGORITHM_H
#define MYALGORITHM_H

#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>
//#include <parallel/algorithm> __gnu_parallel::sort instead of std::sort

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
    void stablesort(std::vector<int>& tosort, std::vector<int>& reorderingvector); /// REMOVE AND USE THE ONE BELOW?
    
    // 'stablesortparallel' sorts in parallel according to tosort[0] then 
    // according to tosrt[1],... ALL int* arrays have size 'numentries'.
    // The output is a reordering vector.
    int* stablesortparallel(std::vector<int*> tosort, int numentries);
    
    std::vector<int> intersect(std::vector<int> a, std::vector<int> b);
};

#endif
