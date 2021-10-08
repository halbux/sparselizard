// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "myalgorithm.h"
#include "polynomials.h"

// Element type numbers are as follows:
//
// 0  --> Point
// 1  --> Line        order 1
// 2  --> Triangle    order 1
// 3  --> Quadrangle  order 1
// 4  --> Tetrahedron order 1
// 5  --> Hexahedron  order 1
// 6  --> Prism       order 1
// 7  --> Pyramid     order 1
// 8  --> Line        order 2
// 9  --> Triangle    order 2
// 10 --> Quadrangle  order 2
// 11 --> Tetrahedron order 2
// 12 --> Hexahedron  order 2
// 13 --> Prism       order 2
// 14 --> Pyramid     order 2
// 15 --> Line        order 3
// ...
//
// Only complete order curved elements are defined.
//
//
// The node, edge and surface ordering in the reference elements is the same as for the GMSH software and is described below.
// The high order nodes for the geometrical element curvature are numbered in the same way as in GMSH.
//
//
// - POINT (element number is 0)
//
// - LINE (element number is 1)
//
// 0----------1 --> ki
//
// (-1) --> Coordinate of 0
// ( 1) --> Coordinate of 1
//
// (0, 1) --> 1st edge
//
// - TRIANGLE (element number is 2)
//
// eta
// ^
// |
// 2
// |`\
// |  `\
// |    `\
// |      `\
// |        `\
// 0----------1 --> ki
//
// (0, 0) --> Coordinate of 0
// (1, 0) --> Coordinate of 1
// (0, 1) --> Coordinate of 2
//
// (0, 1) --> 1st edge
// (1, 2) --> 2nd edge
// (2, 0) --> 3rd edge
//
// (0, 1, 2) --> 1st surface
// 
// - QUADRANGLE (element number is 3)
//
//      eta
//       ^
//       |
// 3-----------2
// |     |     |
// |     |     |
// |     +---- | --> ki
// |           |
// |           |
// 0-----------1
//
// (-1, -1) --> Coordinate of 0
// ( 1, -1) --> Coordinate of 1
// ( 1,  1) --> Coordinate of 2
// (-1,  1) --> Coordinate of 3
//
// (0, 1) --> 1st edge
// (1, 2) --> 2nd edge
// (2, 3) --> 3rd edge
// (3, 0) --> 4th edge
//
// (0, 1, 2, 3) --> 1st surface
//
// - TETRAHEDRON (element number is 4)
//
//                   eta
//                  .
//                ,/
//               /
//            2
//          ,/|`\
//        ,/  |  `\
//      ,/    '.   `\
//    ,/       |     `\
//  ,/         |       `\
// 0-----------'.--------1 --> ki
//  `\.         |      ,/
//     `\.      |    ,/
//        `\.   '. ,/
//           `\. |/   
//              `3
//                 `\.
//                    ` phi
//
// ( 0,  0,  0) --> Coordinate of 0
// ( 1,  0,  0) --> Coordinate of 1
// ( 0,  1,  0) --> Coordinate of 2
// ( 0,  0,  1) --> Coordinate of 3
//
// (0, 1) --> 1st edge
// (1, 2) --> 2nd edge
// (2, 0) --> 3rd edge
// (3, 0) --> 4th edge
// (3, 2) --> 5th edge
// (3, 1) --> 6th edge
//
// (0, 2, 1) --> 1st surface
// (0, 1, 3) --> 2nd surface
// (0, 3, 2) --> 3rd surface
// (3, 1, 2) --> 4th surface
//
// - HEXAHEDRON (element number is 5)
//
//       eta
// 3----------2
// |\     ^   |\
// | \    |   | \
// |  \   |   |  \
// |   7------+---6
// |   |  +-- |-- | -> ki
// 0---+---\--1   |
//  \  |    \  \  |
//   \ |     \  \ |
//    \|     phi \|
//     4----------5
//
// (-1, -1, -1) --> Coordinate of 0
// ( 1, -1, -1) --> Coordinate of 1
// ( 1,  1, -1) --> Coordinate of 2
// (-1,  1, -1) --> Coordinate of 3
// (-1, -1,  1) --> Coordinate of 4
// ( 1, -1,  1) --> Coordinate of 5
// ( 1,  1,  1) --> Coordinate of 6
// (-1,  1,  1) --> Coordinate of 7
//
// (0, 1) --> 1st edge
// (0, 3) --> 2nd edge
// (0, 4) --> 3rd edge
// (1, 2) --> 4th edge
// (1, 5) --> 5th edge
// (2, 3) --> 6th edge
// (2, 6) --> 7th edge
// (3, 7) --> 8th edge
// (4, 5) --> 9th edge
// (4, 7) --> 10th edge
// (5, 6) --> 11th edge
// (6, 7) --> 12th edge
//
// (0, 3, 2, 1) --> 1st surface
// (0, 1, 5, 4) --> 2nd surface
// (0, 4, 7, 3) --> 3rd surface
// (1, 2, 6, 5) --> 4th surface
// (2, 3, 7, 6) --> 5rd surface
// (4, 5, 6, 7) --> 6th surface
//
// - PRISM (element number is 6)
//
//           phi
//            ^
//            |
//            3
//          ,/|`\
//        ,/  |  `\
//      ,/    |    `\
//     5------+------4
//     |      |      |
//     |    ,/|`\    |
//     |  ,/  |  `\  |
//     |,/    |    `\|
//    ,|      |      |\
//  ,/ |      0      | `\
// eta |    ,/ `\    |   ki
//     |  ,/     `\  |
//     |,/         `\|
//     2-------------1
//
// ( 0,  0, -1) --> Coordinate of 0
// ( 1,  0, -1) --> Coordinate of 1
// ( 0,  1, -1) --> Coordinate of 2
// ( 0,  0,  1) --> Coordinate of 3
// ( 1,  0,  1) --> Coordinate of 4
// ( 0,  1,  1) --> Coordinate of 5
//
// (0, 1) --> 1st edge
// (0, 2) --> 2nd edge
// (0, 3) --> 3rd edge
// (1, 2) --> 4th edge
// (1, 4) --> 5th edge
// (2, 5) --> 6th edge
// (3, 4) --> 7th edge
// (3, 5) --> 8th edge
// (4, 5) --> 9th edge
//
// (0, 2, 1)    --> 1st surface
// (3, 4, 5)    --> 2nd surface
// (0, 1, 4, 3) --> 3rd surface
// (0, 3, 5, 2) --> 4th surface
// (1, 2, 5, 4) --> 5rd surface
//
// - PYRAMID (element number is 7)
//
//                4
//              ,/|\
//            ,/ .'|\
//          ,/   | | \
//        ,/    .' | `.
//      ,/      |  '.  \
//    ,/       .'phi|   \
//  ,/         |  ^ |    \
// 0----------.'--|-3    `.
//  `\        |   |  `\    \
//    `\     .'   +----`\ - \ -> eta
//      `\   |    `\     `\  \
//        `\.'      `\     `\`
//           1----------------2
//                     `\
//                       ki
//
// (-1, -1,  0) --> Coordinate of 0
// ( 1, -1,  0) --> Coordinate of 1
// ( 1,  1,  0) --> Coordinate of 2
// (-1,  1,  0) --> Coordinate of 3
// ( 0,  0,  1) --> Coordinate of 4
//
// (0, 1) --> 1st edge
// (0, 3) --> 2nd edge
// (0, 4) --> 3rd edge
// (1, 2) --> 4th edge
// (1, 4) --> 5th edge
// (2, 3) --> 6th edge
// (2, 4) --> 7th edge
// (3, 4) --> 8th edge
//
// (0, 1, 4)    --> 1st surface
// (3, 0, 4)    --> 2nd surface
// (1, 2, 4)    --> 3rd surface
// (2, 3, 4)    --> 4th surface
// (0, 3, 2, 1) --> 5th surface
//
//
// CREDITS:
//
// ASCII ART 3D element representation comes from the open-source GMSH meshing tool documentation.
//
// Unless explicitly mentionned with 'curved'
// in the name, everything refers to the straight 
// version of the actual curved element 

class element
{

    private:
    
        int curvedtypenumber = -1;
        std::vector<int> curvednodelist = {};
        
        polynomials mypolynomials;
        
        // 'getnodesinsurface' does the actual work of 'getnodesintriangle' (below) if the
        // second argument is true and the third is false and of 'getnodesinquadrangle'
        // if the second is false and the third is true. No error check is performed here.
        std::vector<int> getnodesinsurface(int surfaceindex, bool faceistriangle, bool faceisquadrangle);
        
        std::vector<std::vector<int>> splitline(int splitnum);
        std::vector<std::vector<int>> splittriangle(int splitnum, std::vector<int>& edgenumbers);
        std::vector<std::vector<int>> splitquadrangle(int splitnum);
        std::vector<std::vector<int>> splittetrahedron(int splitnum, std::vector<int>& edgenumbers);

        // Get the node-connectivity of a split TRIANGLE element.
        void getsplitconnectivity(std::vector<bool>& connectivity, std::vector<std::vector<int>>& splitdefinition);
        // Get the node-connectivity of a TETRAHEDRON based on the node-connectivity of its faces:
        void getsplitconnectivity(std::vector<bool>& volumeconnectivity, std::vector<std::vector<bool>>& faceconnectivity);
        
        void deducetets(std::vector<bool>& connectivity, std::vector<std::vector<bool>>& tetdefs, int originnode, int numinloop, int curnode, std::vector<bool> isnodeused);
        void deducetets(std::vector<bool>& connectivity, std::vector<std::vector<bool>>& tetdefs);
        
    public:
    
        element(void) {};
        // Set the element name:
        element(std::string elementname);
        // Set the curved type number:
        element(int number);
        // Set the type number and curvature order:
        element(int number, int curvatureorder);
        
        void setnodes(std::vector<int>& nodelist);
        
        std::vector<int> getnodes(void);
        std::string gettypename(void);                                    
        // 'gettypenameconjugation' is singular for 0 or 1 as input, plural otherwise.
        std::string gettypenameconjugation(int numberofelements);                
        bool iscurved(void);                                        
        int getcurvatureorder(void);                                
        // Get the straight type number corresponding to the curved element:
        int gettypenumber(void);                        
        int getcurvedtypenumber(void);                                    
        int countcurvednodes(void);            
        int getelementdimension(void);    
        
        // Number of elements of type typenum/of dimension dim in the element. The curvature nodes are not counted.
        int counttype(int typenum);
        int countdim(int dim);
                        
        int countnodes(void);                        
        int countedges(void);    
        int countfaces(void);                                        
        int counttriangularfaces(void);                    
        int countquadrangularfaces(void);                
        int countvolumes(void);
        
        // Returns true if point with coordinates (ki, eta, phi) is inside the reference element:
        bool isinsideelement(double ki, double eta, double phi);
        // Check if a set of coordinates is inside a straight element (corner nodes must be in 
        // usual element ordering). Only lines in 1D, faces in 2D and volumes in 3D are accepted.
        // Provide an absolute roundoff noise threshold.
        void isinsideelement(std::vector<double>& coords, std::vector<double>& cornercoords, std::vector<bool>& isinside, double roundoffnoise);
        
        // Return the corner node/edge/face index at which each reference coordinate is (return -1 if at none).
        // Can only be called for elements of higher dimension than the requested object.
        void atnode(std::vector<double>& refcoords, std::vector<int>& nodenums, double roundoffnoise = 1e-8);
        void atedge(std::vector<double>& refcoords, std::vector<int>& edgenums, double roundoffnoise = 1e-8);
        void atface(std::vector<double>& refcoords, std::vector<int>& facenums, double roundoffnoise = 1e-8);
        
        // Get the barycenter of each edge/face of a straight element whose node coordinates are provided:
        std::vector<double> getedgebarycenter(std::vector<double>& nc);
        std::vector<double> getfacebarycenter(std::vector<double>& nc);
        
        // Gives the length of the reference line, surface of the reference 
        // triangle/quadrangle and volume of the reference volume elements.
        double measurereferenceelement(void);
        
        // 'istriangularface(i)' returns true if the ith face of the element is triangular
        // and false otherwise (e.g. when quadrangular). 
        bool istriangularface(int facenum);
        // For prisms, to know if the edge is horizontal or vertical:
        bool ishorizontaledge(int edgenum);
        
        // 'getnodesintriangle(i)' gives the ordered node number list 
        // forming the ith triangular face in the element object. 
        // An error is returned if none is found (e.g. for a line element).
        // The node, edge and surface ordering is defined at the top of this header.
        std::vector<int> getnodesinline(int lineindex);
        std::vector<int> getnodesintriangle(int triangleindex);
        std::vector<int> getnodesinquadrangle(int quadrangleindex);
        
        // 'getedgesdefinitionbasedonnodes()[2*i+j]' gives the index of the jth node
        // in the ith edge [node1edgei node2edgei] in the element. The edges ordering 
        // is as described at the top of this file. Only the corner nodes are given,
        // not the curvature-related inner nodes.
        std::vector<int> getedgesdefinitionsbasedonnodes(void);        
        // In 'getfacesdefinitionsbasedonnodes' first all triangular faces 
        // are listed then the quadrangular ones.
        std::vector<int> getfacesdefinitionsbasedonnodes(void);
        std::vector<int> getfacesdefinitionsbasedonedges(void);                                            

        bool iselementedgeorface(void);        

        // This function outputs the coordinates in format {x1,y1,z1,x2,y2,z2,x3,...}
        // of the phi then eta then ki coordinate-ordered nodes in the curved element.
        std::vector<double> listnodecoordinates(void);
        
        // Get the straight element type number corresponding to a number of nodes and dimension:
        int deducetypenumber(int elemdim, int numnodes);
        
        // Calculate the physical coordinates corresponding to a set of reference coordinates.
        // The coordinates in 'nodecoords' are taken from index 'firstindex' and following.
        // If 'returnnodecoords' is true then the node coordinates are returned as they are.
        std::vector<double> calculatecoordinates(std::vector<double>& refcoords, std::vector<double>& nodecoords, int firstindex = 0, bool returnnodecoords = false);
        
        // Count the number of elements of each type in the n times full-split element:
        std::vector<int> fullsplitcount(int n);
        // Full-split n times multiple elements whose node coordinates are provided as argument.
        void fullsplit(int n, std::vector<std::vector<double>>& splitcoords, std::vector<double>& unsplitcoords);
        
        // Get the full-split subelement definition based on their corner reference coordinates.
        // For tetrahedra provide the through-edge to be used in the cut (choose 0 for edge 0-4, 1 for 1-3, 2 for 2-5).
        void fullsplit(std::vector<std::vector<double>>& cornerrefcoords, int throughedgenum = -1);
        
        // Select the through-edge to get the best quality tetrahedron split.
        int choosethroughedge(std::vector<double>& nodecoords);
        
        // Get the reference node numbers defining the transition element split:
        std::vector<std::vector<int>> split(int splitnum, std::vector<int>& edgenumbers);
        // Get the reference coordinates corresponding to the reference node numbers:
        void numstorefcoords(std::vector<int>& nums, std::vector<double>& refcoords);
        
        // Write to file multiple elements with provided coordinates for debug:
        void write(std::string filename, std::vector<double> coords);
        
};

#endif
