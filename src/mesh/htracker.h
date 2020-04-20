// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object manages the mesh h-adaptivity for all element types.

#ifndef HTRACKER_H
#define HTRACKER_H

#include <iostream>
#include <vector>
#include <string>
#include "element.h"

class htracker
{

    private:
        
        // Store in compressed format all split info needed:
        // 
        // 1. 1/0 for fullsplit/transition element (if 0 this element stops here)
        // 3. 00/01/10 for 0/1/2 through-edge tetrahedron fullsplit (only for tetrahedra)
        // 4. 1/0 for fullsplit/transition of first subelement
        // 5. ...
        // 6. 1/0 for fullsplit/transition of second subelement
        // 7. ...
        //
        // Elements types are ordered from lowest type number to highest.
        //
        std::vector<bool> splitdata = {};
        
        int maxdepth = 0;
        int numleaves = 0;
        
        // Number of elements in each type at the 0 depth:
        std::vector<int> originalcount = {};
        
        // Number of subelements in each fullsplit type:
        std::vector<int> numsubelems = {1,2,4,4,8,8,8,10};
        
        // Information needed by the leaf cursor:
        int cursorposition = -1;
        int currentdepth = -1;
        int curtypeorigcountindex = -1;
        std::vector<int> parenttypes = {};
        std::vector<int> indexesinclusters = {};
        
    public:

        // Number of highest dimension elements for each type (provide a vector of length 8).
        htracker(std::vector<int> numelemspertype);
        
        int countleaves(void);
        int getmaxdepth(void);
        
        // Place cursor before first leaf:
        void resetcursor(void);
        // Move cursor forward (crashes when exceeding number of leaves). Position might not be at a leaf.
        // The through-edge number needed for the split is returned (-1 if none or no split).
        int next(void);
        // Check if the cursor is at a leaf:
        bool isatleaf(void);
        // Count the number of splits used to reach the current position:
        int countsplits(void);
        // Get the element type number at current position:
        int gettype(void);
        // Get the element type number of the current parent (-1 if at depth 0):
        int getparenttype(void);
        // Get the index of the current element in the cluster:
        int getindexincluster(void);
        
        // Count the number of splits used to reach each leaf (modifies the cursor):
        void countsplits(std::vector<int>& numsplits);
        // Get the element type number of all leaves (modifies the cursor):
        void gettype(std::vector<int>& types);
        // Count the number of elements of each type after adaptation (modifies the cursor):
        std::vector<int> countintypes(void);
        
        // Group/keep/split (-1/0/1) the requested leaves. Leaves are grouped if at least one
        // is tagged for grouping and no other is split or tagged for splitting in the same cluster.
        // The argument vectors must have a size equal to the number of leaves.
        void adapt(std::vector<int>& operations, std::vector<int>& throughedgenums);
        
        // Print the tree:
        void print(void);
        
        // Get length of 'splitdata' vector:
        int countbits(void);
        
        
        // Get the reference coordinates in the original elements of the adapted elements corner nodes:
        void getadaptedrefcoords(std::vector<std::vector<double>>& arc);
    
};

#endif

