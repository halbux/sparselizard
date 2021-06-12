// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object manages the decomposition of the mesh into domains.

#ifndef DTRACKER_H
#define DTRACKER_H

#include <iostream>
#include <vector>
#include "rawmesh.h"
#include "slmpi.h"


class dtracker
{

    private:
    
        // Skin of the global geometry (-1 if empty):
        int myglobalgeometryskin = -1;
        
        // Number of overlap layers (0 for no-overlap):
        int mynumoverlaplayers = -1;

        std::weak_ptr<rawmesh> myrawmesh;

        // Neighbours (without this rank and without duplicates, sorted ascendingly):
        std::vector<int> myneighbours = {};
        // Direct access (length is numranks):
        std::vector<bool> myisneighbour = {};
    
        // No-overlap interfaces (length is 3*numranks, -1 if none).
        // Entry 3*r+i is the interface of i-dimensional elements with rank r:
        std::vector<int> mynooverlapinterfaces = {};

        // Inner and outer overlap interfaces (length is 3*numranks, -1 if none).
        // Entry 3*r+i is the interface of i-dimensional elements with rank r:
        std::vector<int> myinneroverlapinterfaces = {};
        std::vector<int> myouteroverlapinterfaces = {};
        // Inner and outer overlaps and their skins (length is numranks):
        std::vector<int> myinneroverlaps = {};
        std::vector<int> myouteroverlaps = {};
        std::vector<int> myinneroverlapskins = {};
        std::vector<int> myouteroverlapskins = {};
        
        // Map entry [n][type][i] gives the element number of type 'type' in this domain that corresponds to the
        // ith element of that type in the nth neighbour domain. Only elements in the outer-overlap/no-overlap
        // interfaces can be mapped. Curvature nodes cannot be mapped. The map has value -1 by default.
        std::vector<std::vector<std::vector<int>>> mymaptothisdomain = {};
        
        // Global node numbers for every domain node (unrelated to the local node numbers, -1 for curvature nodes):
        std::vector<long long int> myglobalnodenumbers = {};
        

        // Discover up to 'numtrialelements' neighbours that share cell-1 dimension elements with this rank.
        // The barycenter of all elements shared with the neighbours to discover must be provided as argument.
        // The number of barycenters provided on all ranks is returned (if all zero an empty vector is returned).
        std::vector<int> discoversomeneighbours(int numtrialelements, std::vector<double>& interfaceelembarys, std::vector<int>& neighboursfound);

        // Upon return 'inneighbours[i]' is the number of the neighbour touching the ith interface element (-1 if no neighbour touching).
        // The neighbours provided must be unique and sorted ascendingly. The number of interface elements for each rank must be provided in 'allnumelementsininterface'.
        void discoverinterfaces(std::vector<int> neighbours, std::vector<double>& interfaceelembarys, std::vector<int>& allnumelementsininterface, std::vector<int>& inneighbour);

        // Find new interfaces and populate the output accordingly. Return false if no more interfaces can be found on any rank.
        bool discovercrossinterfaces(std::vector<int>& interfacenodelist, std::vector<int>& interfaceedgelist, std::vector<std::vector<bool>>& isnodeinneighbours, std::vector<std::vector<bool>>& isedgeinneighbours);

        // Define the inner overlaps and their skins:
        void defineinneroverlaps(void);
        // Define the outer overlaps and their skins:
        void exchangeoverlaps(void);
        // Exchange the physical regions on the overlaps:
        void exchangephysicalregions(void);
        // Define the outer overlap interfaces.
        // This works for any geometry whose global skin region does not intersect itself.
        void defineouteroverlapinterfaces(void);
        // Define the inner overlap interfaces:
        void defineinneroverlapinterfaces(void);

        // Map the outer-overlap/no-overlap interfaces:
        void mapnooverlapinterfaces(void);
        void mapoverlapinterfaces(void);
        
        // Create the global node numbers for every domain node:
        void createglobalnodenumbersnooverlap(void);
        void createglobalnodenumbersoverlap(void);

    public:

        void errorundefined(void);

        dtracker(std::shared_ptr<rawmesh> rm, int globalgeometryskin, int numoverlaplayers);

        std::shared_ptr<rawmesh> getrawmesh(void);
        
        bool isoverlap(void);
        
        // Set manually the no-overlap connectivity of this rank.
        // 'nooverlapinterfaces[3*i+j]' is the interface of j-dimensional elements with the ith neighbour (-1 if none).
        void setconnectivity(std::vector<int>& neighbours, std::vector<int>& nooverlapinterfaces);
        // Discover automatically the no-overlap connectivity of this rank.
        // This works for any geometry whose global skin region does not intersect itself.
        void discoverconnectivity(int numtrialelements = 10, int verbosity = 0);

        // Exchange the overlaps (and the physical regions included) with the neighbours and define the overlap interfaces:
        void overlap(void);
    
        // Map the outer-overlap/no-overlap interfaces:
        void mapinterfaces(void);
        
        // Create the global node numbers for every domain node:
        void createglobalnodenumbers(void);
        
        int countneighbours(void);
        std::vector<int> getneighbours(void);
        int getneighbour(int neighbourindex);

        bool isneighbour(int neighbour);
        
        // Return -1 if not defined:
        int getnooverlapinterface(int neighbour, int elementdimension);
        int getinneroverlapinterface(int neighbour, int elementdimension);
        int getouteroverlapinterface(int neighbour, int elementdimension);
        int getinneroverlap(int neighbour);
        int getouteroverlap(int neighbour);
        int getinneroverlapskin(int neighbour);
        int getouteroverlapskin(int neighbour);

        std::vector<std::vector<std::vector<int>>>* getmap(void);
        
        long long int* getglobalnodenumbers(void);
        void writeglobalnodenumbers(std::string filename);

        // Return the list of all ddm-specific regions:
        std::vector<int> listddmregions(void);
        std::vector<bool> isddmdisjointregion(void);
        
        // Know which disjoint region is in the no-overlap region/is owned by this rank:
        std::vector<bool> isdisjointregioninnooverlap(void);
        std::vector<bool> isdisjointregionowned(void);
        
        // Print connectivity information:
        void print(void);
        // Write no-overlap domain interfaces:
        void writeinterfaces(std::string filename);

};

#endif

