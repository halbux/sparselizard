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
#include <memory>
#include "element.h"
#include "rawmesh.h"

class rawmesh;

class htracker
{

    private:
    
        // Original mesh on which this tracker is based:
        std::weak_ptr<rawmesh> myoriginalmesh;
    
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
        // At [i][2*j+0] and [i][2*j+1] this gives the type and index of the corresponding original element:
        std::vector<std::vector<int>> originalsoftransitions = {};
        
        // Transition elements renumbering:
        std::vector<std::vector<int>> touser = {};
        std::vector<std::vector<int>> toht = {};

    public:

        htracker(void) {};
        // Provide the original 'rawmesh' object (or the latter two arguments for debug):
        htracker(std::shared_ptr<rawmesh> origmesh, int curvatureorder = -1, std::vector<int> numelemspertype = {});
        
        std::shared_ptr<rawmesh> getoriginalmesh(void);
        
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
        
        // Get the split depth of all leaves:
        void countsplits(std::vector<int>& numsplits);
        // Get the element type number of all leaves:
        void gettype(std::vector<int>& types);
        // Get the original element number of all leaves:
        void getoriginalelementnumber(std::vector<int>& oen);
        // Get the original element type and index in type of all leaves:
        void getoriginalelement(std::vector<int>& oet, std::vector<int>& oei);
        // For each original element get the number of leaves of each type (by blocks of 8).
        void countsons(std::vector<int>& numsons);
        
        // Count the number of leaves of each type:
        std::vector<int> countintypes(void);
        // Count the number of transition elements of a type:
        int counttransitions(int elementtypenumber);
        
        // Update the operations to how they will actually be treated (priority is split > noop > group):
        void fix(std::vector<int>& operations);
        
        // Group/keep/split (-1/0/1) the requested leaves. Leaves are grouped if all
        // leaves in the cluster are tagged for grouping and none is already split.
        // NEIGBOUR ELEMENTS CANNOT DIFFER BY MORE THAN ONE SPLIT LEVEL.
        // The argument vector must have a size equal to the number of leaves.
        void adapt(std::vector<int>& operations);    
        
        // Get the reference/physical corner node coordinates (physical is optional) of all leaves.
        // Leaves of each type are ordered in the way the appear in the tree.
        void atleaves(std::vector<std::vector<double>>& arc, std::vector<std::vector<double>>& apc, bool withphysicals);
        
        // Get the node coordinates 'ac' of all transition elements after adaptation.
        void getadaptedcoordinates(std::vector<std::vector<double>>& ac);
    
        // Get an upper bound of number of transition elements that can be created from all leaves:
        std::vector<int> countupperbound(void);
    
        // Renumber the transition elements:
        void renumbertransitions(std::vector<std::vector<int>>& renumbering);
        
        // Get the reference coordinate in the original element 'orc' corresponding to the reference coordinates in the transition elements 'rc'.
        // 'ad[t][i]' gives the first position in 'rc' where the ith transition element of type t is. 'oad[i]' does that for the first position
        // of the ith original element in 'orc'. 'ad' and 'oad' have a 1 longer size and their last value is the vec size.
        // 'maprctoorc[i][j]'x3 gives the position in 'orc' corresponding to the point at 'rc[i][3*j]'. 
        // DO NOT CALL THIS YOURSELF (RENUMBERING IS NOT TAKEN INTO ACCOUNT).
        void tooriginal(std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& oad, std::vector<double>& orc, std::vector<std::vector<int>>& maprctoorc);
        // Inverse function of the above. 'maporctorc[2*i+0]' and 'maporctorc[2*i+1]' respectively 
        // give t and j such that the ith point in 'orc' corresponds to 'rc[t][3*j]'. 
        // DO NOT CALL THIS YOURSELF (RENUMBERING IS NOT TAKEN INTO ACCOUNT).
        void fromoriginal(std::vector<int>& oad, std::vector<double>& orc, std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& maporctorc);
        
        // Get the corresponding point at the target. To each rc[i][3*j] is associated its corresponding transition element type and index in type
        // (in targettranselems[i][2*j+0] and targettranselems[i][2*j+1]) and the associated reference coordinate in the target transition element. 
        void getattarget(std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, htracker* target, std::vector<std::vector<int>>& targettranselems, std::vector<std::vector<double>>& targetrefcoords);

        // For a given transition element get the node/edge/face indexes of the original element (respectively in 'on', 'oe', 'of')
        // in which each node/edge/face is (-1 if none). Data is only returned for boundaries.
        void atoriginal(int transitiontype, int transitionnumber, int& originaltype, int& originalnumber, std::vector<int>& on, std::vector<int>& oe, std::vector<int>& of);
    
        // Take a POSITIVE value for each leaf at the origin and return one for each leaf at the target. When values have to
        // be merged the highest value is selected. This and the target htracker cannot differ by more than one adaptation. 
        void getattarget(std::vector<int>& olv, htracker* target, std::vector<int>& tlv);
             
        // Get the leaf number of a transition element:
        int getleafnumber(int transitiontype, int transitionnumber);
         
        // Print the tree:
        void print(void);
        
        // Length of 'splitdata':
        int countbits(void);
};

#endif

