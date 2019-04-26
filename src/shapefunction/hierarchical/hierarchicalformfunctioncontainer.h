// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef HIERARCHICALFORMFUNCTIONCONTAINER_H
#define HIERARCHICALFORMFUNCTIONCONTAINER_H

#include <iostream>
#include <vector>
#include <string>
#include "densematrix.h"
#include "polynomial.h"
#include "orientation.h"

using namespace std;

class hierarchicalformfunctioncontainer
{

	private:
	
        std::string myformfunctiontypename;
        int myelementtypenumber;
        std::vector<double> myevaluationpoints;
        
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
		
	public:
	
        hierarchicalformfunctioncontainer(void) {};
        hierarchicalformfunctioncontainer(std::string formfunctiontypename, int elementtypenumber, std::vector<double> evaluationpoints);
        
        // Know the highest order available in the container.
		int gethighestorder(void) { return val.size()-1; };
		
		// 'set' adds the evaluated form function polynomial to the val container. 
        // All other int arguments are the same as detailed for 'val' above.
        // This function automatically preallocates the 'val' container.
		void set(int h, int i, int j, int k, int l, int n, polynomial& poly);
		
		// 'tomatrix' puts all form function values corresponding to the 
        // input arguments into a 'densematrix' object. The columns of the 
        // dense matrix correspond to the evaluation points while the rows 
        // correspond to all form functions. The form functions are ordered 
        // in the way defined in 'hierarchicalformfunctioniterator'.
		densematrix tomatrix(int totalorientation, int order, int whichderivative, int component);
        
        // Print all form function values for debug.
        void print(bool printallderivatives);
		
};

#endif
