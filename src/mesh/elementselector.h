// sparselizard - Copyright (C) 2017- A. Halbach, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef ELEMENTSELECTOR_H
#define ELEMENTSELECTOR_H

#include <iostream>
#include "disjointregions.h"
#include "myalgorithm.h"
#include "elements.h"
#include "universe.h"
#include <vector>

class elementselector
{
	private:
        
        elementselector(void) {};
	
        // The disjoint regions from which the elements originate.
        std::vector<int> mydisjointregionnumbers;
	
        // The current selected total orientation:
        int currenttotalorientation;
        
        // The index range in 'elems' of the current total orientation.
		int currentrangebegin = 0;
        int currentrangeend = 0;
        
        // The currently selected disjoint regions.
        std::vector<int> selecteddisjointregions = {};
        
        // The three containers below are sorted first with increasing 
        // orientation number (if orientation matters) then according to 
        // their disjoint regions in the 'mydisjointregionnumbers' order.
        
        // All element numbers:
		std::vector<int> elems;
        // Their total orientation:
        std::vector<int> totalorientations;
        // Their disjoint region:
        std::vector<int> disjointregions;
        
	public:
	
        // All disjoint region numbers must correspond to a same element type.
		elementselector(std::vector<int> disjointregionnumbers, bool isorientationdependent = true);
		// Select a single element:
        elementselector(int disjointregionnumber, int elemindex);
        
        int getelementdimension(void);
        int getelementtypenumber(void);
        std::vector<int> getdisjointregions(void) { return mydisjointregionnumbers; };
        int gettotalorientation(void) { return currenttotalorientation; };
        
        // Select the next total orientation. Returns false if there is none.
		bool next(void);
        
        // Count number of elements in current orientation:
        int countincurrentorientation(void) { return currentrangeend - currentrangebegin + 1; };
        
        
        // Select specific disjoint regions. Set to {} for no selection. The 
        // selected disjoint regions must all be in 'mydisjointregionnumbers'. 
        void selectdisjointregions(std::vector<int> disjregs);
        
        // Count the elements that are at the same time in the current 
        // total orientation and in the selected disjoint regions.
        int countinselection(void);
		// Get the numbers of all elements that are at the same time in the 
        // current total orientation and in the selected disjoint regions.
		std::vector<int> getelementnumbers(void);
		// Get the indexes (w.r.t. the current total orientation) of all 
        // elements that are at the same time in the current total 
        // orientation and in the selected disjoint regions.
		std::vector<int> getelementindexes(void);

        // Extract a new element selector from the selection.
        elementselector extractselection(void);
        
};

#endif
