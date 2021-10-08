// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// The hierarchical form functions used change with the edges and faces orientations.
// This can lead to field discontinuity at the element interfaces if the orientations are not
// computed in a coherent way. The coherent way used is described below.
//
// EDGE ORIENTATION:
//
// An edge with nodes numbers [i j] has orientation 
// 
// - 0 if i > j | [0 1]
// - 1 if j > i | [1 0]
// 
// The vector on the right of '|' gives the node reordering to be in orientation 0.
//
// TRIANGULAR SURFACE ORIENTATION:
//
// For the orientation of a triangle there are only 2 things of importance:
// - which node has the biggest number
// - what is the > < relation between the second and the third node
//
// A triangular face with node numbers [i j k] has orientation
//
// - 0 if i > j > k | [0 1 2]
// - 1 if i > k > j | [0 2 1]
// - 2 if j > k > i | [1 2 0]
// - 3 if j > i > k | [1 0 2]
// - 4 if k > i > j | [2 0 1]
// - 5 if k > j > i | [2 1 0]
//
// The vector on the right of '|' gives the node reordering to be in orientation 0.
//
// QUADRANGULAR SURFACE ORIENTATION:
//
// For the orientation of a quadrangle there are only 2 things of importance:
// - which node has the biggest number
// - what is the > < relation between the second and the fourth node
// The third node is fixed since it must be at the corner opposite to the first node.
//
// A quadrangular face with node numbers [i j k l] - a node and its next-next node must be at opposite corners - has orientation
//
// - 0 if max(i,j,k,l) = i and j > l | [0 1 2 3]
// - 1 if max(i,j,k,l) = i and l > j | [0 3 2 1]
// - 2 if max(i,j,k,l) = j and k > i | [1 2 3 0]
// - 3 if max(i,j,k,l) = j and i > k | [1 0 3 2]
// - 4 if max(i,j,k,l) = k and l > j | [2 3 0 1]
// - 5 if max(i,j,k,l) = k and j > l | [2 1 0 3]
// - 6 if max(i,j,k,l) = l and i > k | [3 0 1 2]
// - 7 if max(i,j,k,l) = l and k > i | [3 2 1 0]
//
// The vector on the right of '|' gives the node reordering to be in orientation 0.
//
// TOTAL ORIENTATION:
//
// The total orientation of an element is an integer that uniquely identifies all edges and faces
// orientations in the element. The number is computed as follows:
// totalorientation = sum on all edges(orientation of edge i * 2^i) + 2^numedges * sum on all faces(orientation of face j * 8^j).
// E.g. for a quad with edges orientations [0 1 1 0] and face orientation 7 we have:
// totalorientation = 0*2^0 + 1*2^1 + 1*2^2 + 0*2^3   +   2^4 * 7*8^0;

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include <iostream>
#include <vector>
#include "math.h"
#include "element.h"

namespace orientation
{
    // Only the corner nodes must be provided in the node list.
    std::vector<int> gettotalorientation(int elementtypenumber, std::vector<long long int>& nodelist); 
    
    // Count the number of orientations for a given element type:
    int countorientations(int elemtypenum);
    
    // 'getedgesorientationsfromtotalorientation'returns a vector whose 
    // index i gives the orientation of edge number i.
    std::vector<int> getedgesorientationsfromtotalorientation(int totalorientation, int elementtypenumber);
    std::vector<int> getfacesorientationsfromtotalorientation(int totalorientation, int elementtypenumber);
    
    // 'getreorderingtoreferenceedgeorientation' returns a vector of vectors 'a' 
    // such that if n is a vector containing the node list of an edge of 
    // orientation i then n(a(i)) is an edge of orientation number 0. 
    std::vector<std::vector<int>> getreorderingtoreferenceedgeorientation(void);
    std::vector<std::vector<int>> getreorderingtoreferencetriangularfaceorientation(void);
    std::vector<std::vector<int>> getreorderingtoreferencequadrangularfaceorientation(void);
    
    // 'getedgesorientationsinelement' outputs a vector listing the 
    // orientation of all edges in the element.
    std::vector<int> getedgesorientationsinelement(int elementtypenumber, std::vector<long long int>& nodelist);
    std::vector<int> getfacesorientationsinelement(int elementtypenumber, std::vector<long long int>& nodelist);
    
    // 'getorientationofedgeargument' gives the orientation 
    // of the edge whose nodes are provided as input argument
    int getorientationofedge(std::vector<long long int>& physicalnodesinedge);
    int getorientationoftriangle(std::vector<long long int>& physicalnodesintriangle);
    int getorientationofquadrangle(std::vector<long long int>& physicalnodesinquadrangle);
};

#endif
