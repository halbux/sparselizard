// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object is a wrapper of the actual shape object 'rawshape' pointed
// by 'rawshapeptr' and to which most of the functions are redirected.
// The purpose of this object is to wrap the 'rawshape' object for a 
// convenient user experience.
// An additional advantage is that the std::shared_ptr type pointer ensures
// the pointed 'rawshape' object is always available when there is
// at least one shape using it.


#ifndef SHAPE_H
#define SHAPE_H

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "expression.h"
#include "rawshape.h"

class expression;
class rawshape;

class shape
{   

    private:
        
        std::shared_ptr<rawshape> rawshapeptr = NULL;

        // Give an error on a NULL rawshape pointer:
        void errornullpointer(void);

    public:

        shape(void);
        // Define a shape based on a rawshape pointer:
        shape(std::shared_ptr<rawshape> inputptr);

        // Constructor with the coordinates of all nodes in the shape provided as input (for points and lines only):
        shape(std::string shapename, int physreg, std::vector<double> coords);

        // Constructor based on the coordinates of the corner nodes in the shape:
        shape(std::string shapename, int physreg, std::vector<double> coords, int nummeshpts);
        shape(std::string shapename, int physreg, std::vector<double> coords, std::vector<int> nummeshpts);

        // Constructor based on sub-shapes (not for point shapes):
        shape(std::string shapename, int physreg, std::vector<shape> subshapes, int nummeshpts);
        shape(std::string shapename, int physreg, std::vector<shape> subshapes, std::vector<int> nummeshpts);
        shape(std::string shapename, int physreg, std::vector<shape> subshapes);

        // Define a disk:
        shape(std::string shapename, int physreg, std::vector<double> centercoords, double radius, int nummeshpts);
        shape(std::string shapename, int physreg, shape centerpoint, double radius, int nummeshpts);


        void setphysicalregion(int physreg);

        // Move the shape in the x, y and z direction by a value given in the expression.
        // Only the x, y and z coordinate field are allowed in the expression.
        void move(expression u);

        // Shift the shape in the x, y and z direction.
        void shift(double shiftx, double shifty, double shiftz);
        // Scale the shape in the x, y and z direction.
        void scale(double scalex, double scaley, double scalez);
        // Rotate the mesh by alphax, alphay and alphaz degrees around the x, y and z axis respectively:
        void rotate(double alphax, double alphay, double alphaz);

        shape extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection = {0,0,1});
        std::vector<shape> extrude(std::vector<int> physreg, std::vector<double> height, std::vector<int> numlayers, std::vector<double> extrudedirection = {0,0,1});

        // Duplicate the shape and all subshapes. 
        shape duplicate(void);

        // Return the dimension of the shape (0D, 1D, 2D or 3D):
        int getdimension(void);
        
        // Return the coordinates of the shape mesh nodes:
        std::vector<double> getcoords(void);

        std::string getname(void);

        std::vector<shape> getsons(void);

        int getphysicalregion(void);

         std::shared_ptr<rawshape> getpointer(void);
};



#endif
