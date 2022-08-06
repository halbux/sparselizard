// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <iostream>
#include <vector>
#include "element.h"
#include "nodes.h"
#include "physicalregions.h"
#include "disjointregions.h"
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "orientation.h"
#include "gentools.h"
#include "ptracker.h"

class nodes;

class elements
{

    private:
        
        nodes* mynodes;
        physicalregions* myphysicalregions;
        disjointregions* mydisjointregions;
        
        // The curvature order of all elements. There can only 
        // be a single curvature order in the mesh.
        int mycurvatureorder = -1;
        
        // subelementsinelements[typenum][subtypenum] gives for element of 
        // uncurved type number 'typenum' the list of all subelements of
        // uncurved type 'subtypenum' (0, 1, 2 or 3) it is made of. 
        // [typenum][0] gives the node list in format [nodesinelement1 
        // nodesinelement2 ...], [typenum][1] the line list, [typenum][2] 
        // the triangle list and [typenum][3] the quadrangle list.
        // The element itself is not included, e.g. [0][0] is empty.
        // For nodes all curved nodes are provided. 
        std::vector<std::vector<std::vector<int>>> subelementsinelements = std::vector<std::vector<std::vector<int>>>(8, std::vector<std::vector<int>>(4,std::vector<int>(0)));
        
        // For speedup: number of subelements (nodes, lines, triangles 
        // and quadrangles) in every element.
        std::vector<std::vector<int>> numberofsubelementsineveryelement = std::vector<std::vector<int>>(8,std::vector<int>(4,0));
        
        // indisjointregion[typenum][i] gives the disjoint region number 
        // in which the ith element of type 'typenum' is.
        std::vector<std::vector<int>> indisjointregion = std::vector<std::vector<int>>(8, std::vector<int>(0));
        
        // totalorientations[typenum][i] gives the total orientation 
        // number for the ith element of type 'typenum'.
        // totalorientations[typenum] is empty if not applicable.
        std::vector<std::vector<int>> totalorientations = std::vector<std::vector<int>>(8, std::vector<int>(0));
        
        // barycenters[typenum][3*i] stores the x, y and z coordinates 
        // of the barycenter for the ith element of type 'typenum'.
        // The x, y and z coordinates are concatenated one after the other.
        std::vector<std::vector<double>> barycenters = std::vector<std::vector<double>>(8, std::vector<double>(0));
        
        // sphereradius[typenum][i] gives the smallest radius of the sphere (centered on 
        // the barycenter) that surrounds all nodes of the ith element of type 'typenum'.
        // All nodes of the curved element are considered (but the barycenter is the one of the straight element).
        std::vector<std::vector<double>> sphereradius = std::vector<std::vector<double>>(8, std::vector<double>(0));
        // boxdimensions[typenum][3*i] gives the x, y and z distance between the element barycenter and the
        // smallest box (centered on the barycenter) that surrounds all nodes of the ith element of type 'typenum'.
        // All nodes of the curved element are considered (but the barycenter is the one of the straight element).
        std::vector<std::vector<double>> boxdimensions = std::vector<std::vector<double>>(8, std::vector<double>(0));
        
        // Entries in 'edgesatnodes' from index 'adressedgesatnodes[i]' to 'adressedgesatnodes[i+1]-1' are all 
        // edges touching node i. Only the edge corner nodes (not the curvature nodes) have touching edges.
        std::vector<int> adressedgesatnodes = {};
        std::vector<int> edgesatnodes = {};
        // Create the two vectors above:
        void populateedgesatnodes(void);
        
        // Entries in 'cellsattype[0/1/2/3]' from index 'adresscellsattype[][i]' to 'adresscellsattype[][i+1]-1' are all cells 
        // (highest dimension elements) touching node/edge/tri/quad i and their type number in format {type0,cell0,type1,...}.
        std::vector<std::vector<int>> adresscellsattype = std::vector<std::vector<int>>(4, std::vector<int>(0));
        std::vector<std::vector<int>> cellsattype = std::vector<std::vector<int>>(4, std::vector<int>(0));
        // Create the two vectors above:
        void populatecellsattype(int subtype, std::vector<int>& act, std::vector<int>& ct);

    public:
        
        elements(void) {};
        elements(nodes&, physicalregions&, disjointregions&);
        
        nodes* getnodes(void);
        physicalregions* getphysicalregions(void);
        disjointregions* getdisjointregions(void);
        
        // Add an element defined by its element type number, curvature order and 
        // curved nodes. Return the created element number. Only 'pointsinelements' 
        // is changed. 'elementtypenumber' is the UNCURVED element type number.
        int add(int elementtypenumber, int curvatureorder, std::vector<int>& nodelist);
        
        void cleancoordinatedependentcontainers(void);
        
        // 'getsubelement' returns the number of the 'subelementindex'th 
        // subelement of type 'subelementtypenumber' in element number 
        // 'elementnumber' of type 'elementtypenumber'.
        // E.g. getsubelement(0, 3, 245, 1) returns the second node in quadrangle number 245.
        int getsubelement(int subelementtypenumber, int elementtypenumber, int elementnumber, int subelementindex);
        
        int getdisjointregion(int elementtypenumber, int elementnumber, bool errorifnegative = true);
        int gettotalorientation(int elementtypenumber, int elementnumber);
        
        // Check if a subelement has the same orientation as in its parent element:
        std::vector<bool> isflipped(int subelementtypenumber, std::vector<int>& subelementnumbers, int elementtypenumber, std::vector<int>& elementnumbers);
        
        // 'isinelementlists[i]' is true if the ith element of type 'elementtypenumber' is a (sub-)element of any element list provided.
        // The number of true entries in 'isinelementlists' is returned. NULL element list pointers are ignored.
        int istypeinelementlists(int elementtypenumber, std::vector<std::vector<std::vector<int>>*> elementlists, std::vector<bool>& isinelementlists, bool considercurvaturenodes);
        int istypeindisjointregions(int elementtypenumber, std::vector<bool> isdisjregselected, std::vector<bool>& isinelementlists, bool considercurvaturenodes);
        
        // Return the number of elements of a given type:
        int count(int elementtypenumber);
        // Same as above but returns 0 if the element type is not a cell:
        int countcells(int elementtypenumber);
        // Return the number of elements of a given dimension:
        int countindim(int dim);
        // Count in each type:
        std::vector<int> count(void);
        // Return the curvature order:
        int getcurvatureorder(void);
        
        // Get a vector containing all edges touching node 'nodenumber'.
        // For curvature nodes an empty vector is returned. 
        int countedgesonnode(int nodenumber);
        std::vector<int> getedgesonnode(int nodenumber);
        // Get a vector containing all cells (highest dimension elements)
        // touching node/edge/tri/quad 'subnumber'. Format is {type0,cell0,type1,...}.
        int countcellsontype(int subtype, int subnumber);
        std::vector<int> getcellsontype(int subtype, int subnumber);
        
        // Get the x, y or z coordinate of all nodes in the element 
        // (for xyz respectively set to 0, 1 or 2).
        std::vector<double> getnodecoordinates(int elementtypenumber, int elementnumber, int xyz);
        // Get all coordinates at once:
        std::vector<double> getnodecoordinates(int elementtypenumber, int elementnumber);
        
        // Get the elements in format {type0,elemnum0,type1,...} and reference coordinates at the target disjoint regions.
        // All elements have same reference coordinates at the origin. One reference coordinate per element at target.
        // All target disjoint regions must have elements of same dimension.
        void getrefcoordsondisjregs(int origintype, std::vector<int>& elems, std::vector<double>& refcoords, std::vector<int> targetdisjregs, std::vector<int>& targetelems, std::vector<double>& targetrefcoords);
        
        // Get a pointer to the barycenters[elementtypenumber] vector.
        // The 'barycenters' container is populated for the element type if empty. 
        std::vector<double>* getbarycenters(int elementtypenumber);
        // Get a pointer to the sphereradius[elementtypenumber] vector.
        // The 'sphereradius' container is populated for the element type if empty. 
        std::vector<double>* getsphereradius(int elementtypenumber);
        // Get a pointer to the boxdimensions[elementtypenumber] vector.
        // The 'boxdimensions' container is populated for the element type if empty. 
        std::vector<double>* getboxdimensions(int elementtypenumber);
        
        // Get the barycenter of all elements in the flattened element list:
        void getbarycenters(std::vector<std::vector<int>>* elementlist, std::vector<double>& barycenters);
        // Get the barycenter of all requested elements:
        void getbarycenters(int elementtypenumber, std::vector<int>& elementlist, std::vector<double>& barycenters);
        void getbarycenters(int elementtypenumber, std::vector<int>& elementlist, double* barycenters);
        
        // Get the normal (not normed) to a straight face element:
        std::vector<double> getnormal(int elementtypenumber, int elementnumber);
        
        // Get the highest element dimension available:
        int getdimension(void);
        
        // Count the number of interface element types (1/2/4 for a 1D/2D/3D mesh):
        int countinterfaceelementtypes(void);
        
        // Print elements data for debug:
        void printnumber(void);
        void printsubelements(void);
        void printtotalorientations(void);
        // Visualize elements for debug:
        void write(std::string filename, int elementtypenumber, std::vector<int> elementnumbers, std::vector<int> elementvalues);
        // Visualize edge direction for debug:
        void writeedgedirection(int physreg, std::string filename);
        
        // 'computebarycenters' computes the barycenter coordinates of all 
        // elements of type 'elementtypenumber'. The output vector has format 
        // [coordxelem1 coordyelem1 coordzelem1 coordxelem2 ...].
        // Only the corner element nodes are used.
        std::vector<double> computebarycenters(int elementtypenumber);
        
        // 'removeduplicates' removes the duplicated elements.
        std::vector<int> removeduplicates(int elementtypenumber);

        // 'renumber' updates the element numbers in 'linesinelements', 
        // 'trianglesinelements' and 'quadranglesinelements'. 'renumberingvector' 
        // is such that linesinelementsrenumbered = renumberingvector(linesinelements).
        void renumber(int elementtypenumber, std::vector<int>& renumberingvector);    
        // 'reorder' is similar to 'renumber' except that it reorders the elements
        // using 'reorderingvector' such that orderedelements = elements(reorderingvector,:);
        void reorder(int elementtypenumber, std::vector<int>& reorderingvector);    

        // 'explode' extract the subelements from the existing ones. The 
        // elements created are lines (the edges) and triangles/quadrangles 
        // (the faces) that are part of the original elements. This function 
        // creates duplicated elements that have to be uniqued afterwards.
        void explode(void);
        
        // Follow the elements of highest dimension in 'elementlist' in the order provided.
        // For each element follow the subelements of type 'subtype' in the element-subs order.
        // Each subelement is passed only once. If the last argument is provided then only the
        // subelements that are also subelements of any last argument element list are returned.
        // NULL pointers in the last argument are empty lists. Curvature nodes are NOT returned.
        void follow(std::vector<std::vector<int>>* elementlist, int subtype, std::vector<int>& sublist, std::vector<std::vector<std::vector<int>>*> mustbeinelementlists = {});         
        
        
        // To call only after all renumbering and reordering steps:
        void definedisjointregions(void);
        // Reorder and renumber the elements by disjoint regions.
        void reorderbydisjointregions(void);
        // Same but here the renumbering used is provided in 'elementrenumbering' upon return.
        void reorderbydisjointregions(std::vector<std::vector<int>>& elementrenumbering);
        void definedisjointregionsranges(void);
        
        // Get a vector whose index i is true if node i is a corner node:
        std::vector<bool> iscornernode(void);

        // 'orient' defines 'totalorientations'. Note: 'totalorientations'
        // is untouched in all renumbering and reordering steps and
        // should thus be called last, after all other steps.
        void orient(long long int* noderenumbering = NULL);
        
        
        // Make a full copy of this object (linking objects used are the arguments):
        elements copy(nodes* nds, physicalregions* prs, disjointregions* drs);
        
        // Bring this object in the state corresponding to the target ptracker:
        void toptracker(std::shared_ptr<ptracker> originpt, std::shared_ptr<ptracker> targetpt);
        
        // Merge with this object. The renumbering into the merged object of each
        // (sub)element to merge must be provided (can include duplicates).
        // The 'nodes' object known to this object is modified.
        void merge(elements* elstomerge, std::vector<std::vector<int>>& renumbering, std::vector<int>& numduplicates);
        
        // Merge with this object. New elements are appended. Duplicates in 'elstomerge' are
        // removed. Duplicates in the intersection region between 'elstomerge' and this object
        // are also removed. THIS OBJECT CAN NOT HAVE DUPLICATES on the intersection region.
        // CELL DUPLICATES ARE NOT REMOVED. The 'nodes' object known to this object is modified.
        void merge(std::vector<int> intersectionphysregs, elements* elstomerge);
};

#endif
