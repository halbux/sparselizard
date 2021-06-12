// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef RAWMESH_H
#define RAWMESH_H

#include <string>
#include "nodes.h"
#include "elements.h"
#include "physicalregions.h"
#include "disjointregions.h"
#include "universe.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "wallclock.h"
#include "physicalregion.h"
#include "gmshinterface.h"
#include "nastraninterface.h"
#include "element.h"
#include <memory>
#include "shape.h"
#include "rawshape.h"
#include "regiondefiner.h"
#include "petscmesh.h"
#include "ptracker.h"
#include "htracker.h"
#include "rawfield.h"
#include "dtracker.h"

class dtracker;
class htracker;
class nodes;
class elements;
class shape;

class rawmesh : public std::enable_shared_from_this<rawmesh>
{
    private:
        
        nodes mynodes;
        elements myelements;
        physicalregions myphysicalregions;
        disjointregions mydisjointregions;
        
        int mynumsplitrequested = 0;
        void splitmesh(void);
        regiondefiner myregiondefiner;
        
        int mynumber = 0;
        
        // For domain decomposition:
        std::shared_ptr<dtracker> mydtracker = NULL;
        
        // For p-adaptivity:
        std::shared_ptr<ptracker> myptracker = NULL;
        std::vector<std::tuple<std::weak_ptr<rawfield>, expression, int, int, double>> mypadaptdata = {};
    
        // For h-adaptivity (only for original mesh):
        std::shared_ptr<rawmesh> myhadaptedmesh = NULL;
        std::vector<std::tuple<expression, int, int, double>> myhadaptdata = {}; // only one element or empty if not h-adaptive
        // For the h-adapted mesh:
        std::shared_ptr<htracker> myhtracker = NULL;
        
    public:
        
        // 'readfromfile' hands over to the function reading the format of the mesh file.
        void readfromfile(std::string tool, std::string source);
        // 'writetofile' hands over to the function writing the format of the mesh file.
        void writetofile(std::string);

        // 'removeduplicates' removes the duplicated elements (and nodes).
        void removeduplicates(int lasttypetoprocess = 7);
        
        // 'printcount' prints the number of elements for every type.
        void printcount(void);

        rawmesh(void);
        ~rawmesh(void);
        
        nodes* getnodes(void);
        elements* getelements(void);
        physicalregions* getphysicalregions(void);
        disjointregions* getdisjointregions(void);
        std::shared_ptr<dtracker> getdtracker(void);
        std::shared_ptr<ptracker> getptracker(void);
        std::shared_ptr<htracker> gethtracker(void);
        int getmeshnumber(void) { return mynumber; };
        
        // Get a full copy of this rawmesh:
        std::shared_ptr<rawmesh> copy(void);
        
        // Get a full copy of this rawmesh adapted to the target ptracker.
        // If the target ptracker is identical to 'myptracker' then this object is returned.
        std::shared_ptr<rawmesh> getattarget(std::shared_ptr<ptracker> targetpt);

        // Load from file name:
        void load(std::string name, int globalgeometryskin, int numoverlaplayers, int verbosity);   
        // Load from multiple files:
        void load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity);
        // Load from shape vector:
        void load(std::vector<shape> inputshapes, int globalgeometryskin, int numoverlaplayers, int verbosity);

        // Write to file name:
        void write(std::string name, int verbosity);     
        
        // Split each element in the mesh n times:
        void split(int n);
        
        // Get a bool vector telling if the nodes are in a physical region:
        std::vector<bool> isnodeinphysicalregion(int physreg);
        
        // Move the mesh in the x, y and z direction by a value given in the expression.
        void move(int physreg, expression u);
        // 'shift' translates the mesh in the 'x', 'y' and 'z' direction.
        void shift(int physreg, double x, double y, double z);
        // 'rotate' rotates the mesh by ax, ay and az degrees around the x, y and z axis respectively.
        void rotate(int physreg, double ax, double ay, double az);
        // 'scale' scales the mesh in the 'x', 'y' and 'z' direction.
        void scale(int physreg, double x, double y, double z);
        
        // 'getmeshdimension' gives n for a mesh whose highest element dimension is n.
        int getmeshdimension(void);
        
        // Get the physical regions of a given dimension (use -1 for all).
        std::vector<int> getphysicalregionnumbers(int dim);
        
        // Additional region selection tools. Will become effective only after loading the mesh. Can reuse previous selections!
        void selectskin(int newphysreg, int physregtoskin);
        void selectbox(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit);
        void selectsphere(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius);
        void selectlayer(int newphysreg, int physregtoselectfrom, int physregtostartgrowth, int numlayers);
        void selectexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude);
        void selectanynode(int newphysreg, int physregtoselectfrom);
        

        bool adapthp(int verbosity);
        // Get the POSITIVE values[elementtype][elementnumber] at the target mesh (take the highest value 
        // if values must be merged). This and the target mesh cannot differ by more than one adaptation. 
        void getattarget(std::vector<std::vector<int>>& values, std::shared_ptr<rawmesh> target);
        
        // For p-adaptivity:
        void add(std::shared_ptr<rawfield> inrawfield, expression criterion, int loworder, int highorder, double critrange);
        void remove(rawfield* inrawfield);
        bool adaptp(std::vector<std::vector<std::vector<int>>>& neworders, int verbosity);
        
        // For h-adaptivity:
        bool adapth(std::vector<std::vector<int>>& groupkeepsplit, int verbosity);
        void setadaptivity(expression criterion, int lownumsplits, int highnumsplits, double critrange); // critrange -1 for automatic choice

        // FOR DEBUG. The physical regions are replaced by disjoint regions + 1:
        void writewithdisjointregions(std::string);
        // Print the disjoint region list in every physical region:
        void printphysicalregions(void);
        // Print the physical region list in every disjoint region:
        void printdisjointregions(void);
        // Print the elements in every physical region:
        void printelementsinphysicalregions(bool isdebug = false);
        
        // Check for disconnected lower dimension geometrical entities (e.g. disconnected faces in 3D):
        void errorondisconnecteddisjointregion(void);
        
        std::shared_ptr<rawmesh> gethadaptedpointer(void);
        std::shared_ptr<rawmesh> getoriginalmeshpointer(void);
};


#endif
 
