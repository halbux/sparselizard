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
#include "elements.h"

class htracker
{

    private:
    
        // Original elements on which this tracker is based:
        elements* myoriginalelements = NULL;
    
        // Store in compressed format all split info needed:
        // 
        // 1. 1/0 for fullsplit/transition element (if 0 this element stops here)
        // 2. 00/01/10/11 for 0/1/2/undefined through-edge tetrahedron fullsplit (only for tetrahedra)
        // 3. 1/0 for fullsplit/transition element of first subelement
        // 4. ...
        // 5. 1/0 for fullsplit/transition element of second subelement
        // 6. ...
        //
        // Transition elements are numbered for each type in the order in which they appear in the tree.
        // A transition element can therefore be uniquely identified by its type and index in that type.
        //
        std::vector<bool> splitdata = {};
        
        int maxdepth = 0;
        int numleaves = 0;
        
        // Curvature order of the original elements:
        int originalcurvatureorder;
        // Number of elements in each type at the 0 depth:
        std::vector<int> originalcount = {};
        
        // Number of subelements in each fullsplit type:
        std::vector<int> numsubelems = {1,2,4,4,8,8,8,10};
        // Number of corner/curvature nodes for each element type:
        std::vector<int> nn, ncn;
        // Straight/curved reference coordinates for each element type:
        std::vector<std::vector<double>> straightrefcoords;
        std::vector<std::vector<double>> curvedrefcoords;
        
        // Straight/curved element object for each element type:
        std::vector<element> myelems;
        std::vector<element> mycurvedelems;
        
        // Information needed by the tree cursor:
        bool isrefcalc = false;
        int cursorposition = -1;
        int currentdepth = -1;
        int origindexintype = -1;
        std::vector<int> parenttypes = {};
        std::vector<int> indexesinclusters = {};
        std::vector<std::vector<std::vector<double>>> parentrefcoords = {};
        // End tree cursor information
        
        // Reference coordinates in the original element/leaf number for each transition element.
        // Transition elements of each type are ordered in the way the appear in the tree.
        std::vector<std::vector<double>> transitionsrefcoords = {};
        std::vector<std::vector<int>> leavesoftransitions = {};

    public:

        // Provide the original 'elements' object (or the latter two arguments for debug):
        htracker(elements* origelems, int curvatureorder = -1, std::vector<int> numelemspertype = {});
        
        int countleaves(void);
        int getmaxdepth(void);
        
        // Place cursor at beginning of tree. Request reference coordinate calculations or not.
        void resetcursor(bool calcrefcoords = false);
        // Move cursor to next node (crashes when exceeding number of leaves). Position might not be at a leaf.
        // The through-edge number needed for the split is returned (-1 if n/a, 3 if undefined). In case it
        // has not been defined yet and the 'resetcursor' argument is true it is calculated and stored.
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
        // Get the current reference coordinates in the original element (only if 'calcrefcoords' is true):
        std::vector<double> getreferencecoordinates(void);
        
        // Count the number of splits used to reach each leaf:
        void countsplits(std::vector<int>& numsplits);
        // Get the element type number of all leaves:
        void gettype(std::vector<int>& types);
        // For each original element 'numsons' gives the number of sons of each type (by blocks of 8).
        void countsons(std::vector<int>& numsons);
        // For each leaf get the original element number:
        void getoriginalelementnumber(std::vector<int>& oen);
        
        // Count the number of leaves of each type:
        std::vector<int> countintypes(void);
        
        // Update the operations to how they will actually be treated (this removes individual grouping requests):
        void fix(std::vector<int>& operations);
        
        // Group/keep/split (-1/0/1) the requested leaves. Leaves are grouped if all
        // leaves in the cluster are tagged for grouping and none is already split.
        // NEIGBOUR ELEMENTS CANNOT DIFFER BY MORE THAN ONE SPLIT LEVEL.
        // The argument vector must have a size equal to the number of leaves.
        void adapt(std::vector<int>& operations);    
        
        // Get the reference/physical corner node coordinates (physical is optional) of all leaves.
        // Leaves of each type are ordered in the way the appear in the tree.
        void atleaves(std::vector<std::vector<double>>& arc, std::vector<std::vector<double>>& apc, bool withphysicals);
        
        // Get the node coordinates 'ac' of all transition elements after adaptation. The parent number of each transition element is 
        // provided in 'leafnums' after execution.
        void getadaptedcoordinates(std::vector<std::vector<double>>& ac, std::vector<std::vector<int>>& leafnums, std::vector<double> noisethreshold);
    
        // Get the reference coordinate in the original element 'orc' corresponding to the reference coordinates in the transition elements 'rc'.
        // 'ad[t][i]' gives the first position in 'rc' where the ith transition element of type t is. 'oad[i]' does that for the first position
        // of the ith original element in 'orc'. 'ad' and 'oad' have a 1 longer size and their last value is the vec size.
        void tooriginal(std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& oad, std::vector<double>& orc);
        // Inverse function of the above:
        void fromoriginal(std::vector<int>& oad, std::vector<double>& orc, std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc);
        // + WRITE ONE FCT DOING BOTH ABOVE (arg is htracker&)
        
        // Reduce size for storage:
        void tostorage(void);
    
        // Print the tree:
        void print(void);
        
        // Length of 'splitdata':
        int countbits(void);
};

#endif

