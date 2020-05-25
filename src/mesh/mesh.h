// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object is a wrapper of the actual mesh object 'rawmesh' pointed
// by 'rawmeshptr' and to which most of the functions are redirected.
// The purpose of this object is to wrap the 'rawmesh' object for a 
// convenient user experience.
// An additional advantage is that the std::shared_ptr type pointer ensures
// the pointed 'rawmesh' object is always available when there is
// at least one object using it.

#ifndef MESH_H
#define MESH_H

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "shape.h"
#include "rawmesh.h"
#include "field.h"
#include "expression.h"

class field;
class expression;
class rawmesh;
class shape;

class mesh
{
    private:
    
        bool isloaded = false;
        void errorifloaded(void);
        
        // The actual mesh:
        std::shared_ptr<rawmesh> rawmeshptr = NULL;
    
    public:

        mesh(void);
        mesh(std::string filename, int verbosity = 1, bool legacyreader = true);
        mesh(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity = 1);
        mesh(std::vector<shape> inputshapes, int verbosity = 1);

        // Load from file name:
        void load(std::string name, int verbosity = 1, bool legacyreader = true);   
        // Load from multiple files:
        void load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity = 1);
        // Load from shape vector:
        void load(std::vector<shape> inputshapes, int verbosity = 1);

        // Write to file name:
        void write(std::string name, int verbosity = 1);     
        
        // H-adaptivity:
        bool adapt(int verbosity = 0);
        void setadaptivity(expression criterion, std::vector<field> triggers, int lownumsplits, int highnumsplits, double thresdown = 0.0, double thresup = 0.0, double mincritrange = 0.0);
        void setadaptivity(expression criterion, std::vector<field> triggers, std::vector<double> thresholds, std::vector<int> numsplits, double thresdown = 0.0, double thresup = 0.0, double mincritrange = 0.0);
        
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
        
        // Set this mesh as the one to use:
        void use(void);
        
        std::shared_ptr<rawmesh> getpointer(void) { return rawmeshptr; };
        
};


#endif
 
