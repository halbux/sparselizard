// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef NODES_H
#define NODES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

class nodes
{

    private:
        
        // Coordinates of every node. Format is [x1 y1 z1 x2 y2 z2 ... ].
        std::vector<double> mycoordinates = {};
        
    public:
        
        nodes(void);
        
        // Set the number of nodes.
        void setnumber(int numberofnodes);
        // Get the number of nodes:
        int count(void);
        // Get the coordinates:
        std::vector<double>* getcoordinates(void);
        
        // Print node coordinates for debugging:
        void print(void);
        
        // 'removeduplicates' removes the duplicated nodes in 'mycoordinates'.
        std::vector<int> removeduplicates(void);    
        
        // 'reorder' updates the node orders (i.e. renumbers them) in 
        // 'mycoordinates' based on the input vector. The reordering must be bijective.
        void reorder(std::vector<int>& nodereordering);
        
        // 'getgeometrydimension' gives the max length of the geometry in 
        // the x, y and z dimension for input 0, 1 and 2 respectively.
        double getgeometrydimension(int coord);    
        // Get a round off noise threshold for each coordinate.
        std::vector<double> getnoisethreshold(void);
        
        // Set to zero all negative x coordinates if within the noise range. Throw an error if above the noise range.
        void fixifaxisymmetric(void);
        
};

#endif
