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
        void errorifnotloaded(void);
        
        // The actual mesh:
        std::shared_ptr<rawmesh> rawmeshptr = NULL;
    
    public:

        mesh(void);
        mesh(std::string filename, int verbosity = 1);
        mesh(std::string filename, int globalgeometryskin, int numoverlaplayers, int verbosity = 1);
        mesh(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity = 1);
        mesh(std::vector<shape> inputshapes, int verbosity = 1);
        mesh(std::vector<shape> inputshapes, int globalgeometryskin, int numoverlaplayers, int verbosity = 1);

        // Load from file name:
        void load(std::string name, int verbosity = 1);   
        void load(std::string name, int globalgeometryskin, int numoverlaplayers, int verbosity = 1);   
        // Load from multiple files:
        void load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity = 1);
        // Load from shape vector:
        void load(std::vector<shape> inputshapes, int verbosity = 1);
        void load(std::vector<shape> inputshapes, int globalgeometryskin, int numoverlaplayers, int verbosity = 1);

        // Write to file name:
        void write(std::string name, int verbosity = 1);     
        
        // H-adaptivity:
        void setadaptivity(expression criterion, int lownumsplits, int highnumsplits);
        
        // Split each element in the mesh n times:
        void split(int n = 1);
        
        // Move the mesh in the x, y and z direction by a value given in the expression.
        void move(int physreg, expression u);
        void move(expression u);
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
        void selectskin(int newphysreg, int physregtoskin);
        void selectskin(int newphysreg);
        void selectbox(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit);
        void selectbox(int newphysreg, int selecteddim, std::vector<double> boxlimit);
        void selectsphere(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius);
        void selectsphere(int newphysreg, int selecteddim, std::vector<double> centercoords, double radius);
        void selectlayer(int newphysreg, int physregtoselectfrom, int physregtostartgrowth, int numlayers);
        void selectlayer(int newphysreg, int physregtostartgrowth, int numlayers);
        void selectexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude);
        void selectexclusion(int newphysreg, std::vector<int> physregstoexclude);
        void selectanynode(int newphysreg, int physregtoselectfrom);
        void selectanynode(int newphysreg);
        
        // Set this mesh as the one to use:
        void use(void);
        
        std::shared_ptr<rawmesh> getpointer(void) { return rawmeshptr; };
        
};


#endif
 
