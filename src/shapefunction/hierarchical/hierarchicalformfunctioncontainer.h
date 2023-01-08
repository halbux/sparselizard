// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef HIERARCHICALFORMFUNCTIONCONTAINER_H
#define HIERARCHICALFORMFUNCTIONCONTAINER_H

#include <iostream>
#include <vector>
#include <string>
#include "densemat.h"
#include "polynomial.h"
#include "orientation.h"

using namespace std;

class hierarchicalformfunctioncontainer
{

    private:
    
        std::string myformfunctiontypename;
        int myelementtypenumber;
        
        std::vector<double> myevaluationpoints = {};

        // The values at evaluation points 'myevaluationpoints' of the hierarchical 
        // form functions are stored in the following format:
        //
        // 'val[h][i][j][k][l][m][n][o]' gives the value of the form function
        //
        // - specific to order h (i.e. order 0 is at index 0)
        // - associated to nodes if i is 0, edges if 1, faces if 2 and volume if 3
        // - for node, edge, face or volume number j
        // - for orientation k
        // - number l
        // - for derivative ki if m is 1, eta 2, phi 3, none 0
        // - at component n (0 for x, 1 for y and 2 for z)
        // - evaluation point o
        vector<vector<vector<vector<vector<vector<vector<vector<double>>>>>>>> val = {};
        
        // The form function polynomials are stored in the same format but without [m] and [o]:
        vector<vector<vector<vector<vector<vector<polynomial>>>>>> ffpoly = {};

    public:

        hierarchicalformfunctioncontainer(void) {};
        hierarchicalformfunctioncontainer(std::string formfunctiontypename, int elementtypenumber);
        
        // Know the highest order available in the container.
        int gethighestorder(void) { return val.size()-1; };

        // 'set' adds the form function polynomial to the ffpoly container. 
        // All other int arguments are the same as detailed for 'val' above.
        // This function automatically preallocates the 'val' and 'ffpoly' containers.
        void set(int h, int i, int j, int k, int l, int n, polynomial& poly);
        
        // Evaluate all form function polynomials at the evaluation points provided:
        void evaluate(std::vector<double> evaluationpoints);

        // 'tomatrix' puts all form function values corresponding to the 
        // input arguments into a 'densemat' object. The columns of the 
        // dense matrix correspond to the evaluation points while the rows 
        // correspond to all form functions. The form functions are ordered 
        // in the way defined in 'hierarchicalformfunctioniterator'.
        densemat tomatrix(int totalorientation, int order, int whichderivative, int component);
        densemat tomatrix(int h, int i, int j, int k, int l, int m, int n);

        // Print all form function values for debug.
        void print(bool printallderivatives);

};

#endif
