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
        
        // For p-adaptivity:
        std::shared_ptr<ptracker> myptracker = NULL;
        std::vector<std::tuple<std::weak_ptr<rawfield>, expression,std::vector<double>,std::vector<int>,double,double,double>> mypadaptdata = {};
    
        // For h-adaptivity:
        std::vector<std::shared_ptr<htracker>> myhtracker = {}; // h-tracker to get adapted mesh below
        std::vector<std::shared_ptr<rawmesh>> myhadaptedmesh = {};
        std::vector<std::tuple<expression,std::vector<std::shared_ptr<rawfield>>,std::vector<double>,std::vector<int>,double,double,double>> myhadaptdata = {}; // only one element or empty if not h-adaptive
        std::vector<std::vector<int>> leafnumbersoftransitions = {};
    
    public:
        
        // 'readfromfile' hands over to the function reading the format of the mesh file.
        void readfromfile(std::string);
        // 'writetofile' hands over to the function writing the format of the mesh file.
        void writetofile(std::string);

        // 'sortbybarycenters' sorts the elements (and nodes) according to 
        // the x then y then z coordinates of their barycenter. Round off
        // noise on the coordinates is taken into account in the sorting.
        void sortbybarycenters(void);
        // 'removeduplicates' removes the duplicated elements (and nodes).
        void removeduplicates(void);
        
        // 'printcount' prints the number of elements for every type.
        void printcount(void);

        rawmesh(void);
        
        nodes* getnodes(void);
        elements* getelements(void);
        physicalregions* getphysicalregions(void);
        disjointregions* getdisjointregions(void);
        std::shared_ptr<ptracker> getptracker(void);
        int getmeshnumber(void) { return mynumber; };

        // Load from file name:
        void load(std::string name, int verbosity = 1, bool legacyreader = true);   
        // Load from multiple files:
        void load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity = 1);
        // Load from shape vector:
        void load(std::vector<shape> inputshapes, int verbosity = 1);

        // Write to file name:
        void write(std::string name, int verbosity = 1);     
        
        // Split each element in the mesh n times:
        void split(int n = 1);
        
        // 'shift' translates the mesh in the 'x', 'y' and 'z' direction.
        void shift(int physreg, double x, double y, double z);
        void shift(double x, double y, double z);
        // 'rotate' rotates the mesh by ax, ay and az degrees around the x, y and z axis respectively.
        void rotate(int physreg, double ax, double ay, double az);
        void rotate(double ax, double ay, double az);
        // 'scale' scales the mesh in the 'x', 'y' and 'z' direction.
        void scale(int physreg, double x, double y, double z);
        void scale(double x, double y, double z);
        
        // 'getmeshdimension' gives n for a mesh whose highest element dimension is n.
        int getmeshdimension(void);
        
        // Get the physical regions of a given dimension (use -1 for all).
        std::vector<int> getphysicalregionnumbers(int dim = -1);
        
        // Additional region selection tools. Will become effective only after loading the mesh. Can reuse previous selections!
        void regionskin(int newphysreg, int physregtoskin);
        void boxselection(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit);
        void sphereselection(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius);
        void regionexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude);
        
        // For p-adaptivity:
        void add(std::shared_ptr<rawfield> inrawfield, expression criterion, std::vector<double> thresholds, std::vector<int> orders, double thresdown, double thresup, double mincritrange);
        void remove(rawfield* inrawfield);
        void adaptp(void);
        
        // For h-adaptivity:
        void adapth(void);
        void setadaptivity(expression criterion, std::vector<field> triggers, int lownumsplits, int highnumsplits, double thresdown, double thresup, double mincritrange);
        void setadaptivity(expression criterion, std::vector<field> triggers, std::vector<double> thresholds, std::vector<int> numsplits, double thresdown, double thresup, double mincritrange);

        // FOR DEBUG. The physical regions are replaced by disjoint regions + 1:
        void writewithdisjointregions(std::string);
        // Print the disjoint region list in every physical region:
        void printphysicalregions(void);
        // Print the physical region list in every disjoint region:
        void printdisjointregions(void);
        // Print the elements in every physical region:
        void printelementsinphysicalregions(bool isdebug = false);
        
        std::shared_ptr<rawmesh> getpointer(void);
};


#endif
 
