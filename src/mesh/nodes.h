// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef NODES_H
#define NODES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "myalgorithm.h"

class nodes
{

    private:
        
        // 'roundoffnoiselevel' quantifies the round off noise on the node coordinates:
        double roundoffnoiselevel = 1e-10;
        
        // Coordinates of every node. Format is [x1 y1 z1 x2 y2 z2 ... ].
        std::vector<double> mycoordinates = {};
        
    public:
        
        nodes(void);
        nodes(double roundoffnoise);
        
        // Set the number of nodes.
        void setnumber(int numberofnodes);
        // Get the number of nodes:
        int count(void);
        // Get the coordinates:
        std::vector<double>* getcoordinates(void);
        
        void shift(double xshift, double yshift, double zshift);
        void rotate(double alphax, double alphay, double alphaz);
        
        // Print node coordinates for debugging:
        void print(void);
        
        // 'sortbycoordinates' sorts the nodes according to their x, y and 
        // z coordinates with highest sorting priority for x then y and then z. 
        // 'mycoordinates' is updated accordingly.
        // The output vector v is such that sortedcoordinates(v,:) = coordinates.
        std::vector<int> sortbycoordinates(void);
        // 'removeduplicates' removes the duplicated nodes in 'mycoordinates'.
        // NOTE: 'mycoordinates' must be sorted before the call.
        // The output vector v is such that sortedcoordinates(v,:) = coordinates.
        std::vector<int> removeduplicates(void);    
        
        // 'reorder' updates the node orders (i.e. renumbers them) in 
        // 'mycoordinates' based on the input vector. The input vector is 
        // such that sortedcoordinates = coordinates(nodereordering,:).
        // The reordering must be bijective.
        void reorder(std::vector<int>& nodereordering);
        
        // 'getgeometrydimension' gives the max length of the geometry in 
        // the x, y and z dimension for input 0, 1 and 2 respectively.
        double getgeometrydimension(int coord);    
        // Get a round off noise threshold for each coordinate.
        std::vector<double> getnoisethreshold(void);
        
};

#endif
