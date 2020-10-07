// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object holds a geometrical shape.

#ifndef RAWSHAPE_H
#define RAWSHAPE_H

#include <iostream>
#include <vector>
#include <memory>
#include "mesh.h"
#include "expression.h"

class expression;

class rawshape : public std::enable_shared_from_this<rawshape>
{   

    private:
        
        
    public:

        virtual void move(expression xmove, expression ymove, expression zmove, bool recursively = true);

        virtual void shift(double shiftx, double shifty, double shiftz, bool recursively = true);
        virtual void scale(double scalex, double scaley, double scalez, bool recursively = true);
        virtual void rotate(double alphax, double alphay, double alphaz, bool recursively = true);


        virtual std::shared_ptr<rawshape> extrude(int physreg, double height, int numlayers, std::vector<double> extrudedirection);

        virtual std::shared_ptr<rawshape> duplicate(void);

        // Flip a line direction:
        virtual void flip(void);

        virtual void setphysicalregion(int physreg);

        virtual int getdimension(void);

        virtual std::string getname(void);

        virtual std::vector<std::shared_ptr<rawshape>> getsons(void);

        // Get subshapes (sons are included):
        virtual std::vector<std::shared_ptr<rawshape>> getsubshapes(void);
        // Get all subshapes recursively:
        virtual std::vector<std::shared_ptr<rawshape>> getsubshapesrecursively(void);

        virtual void setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes);
        virtual void setsubshapesrecursively(std::vector<std::shared_ptr<rawshape>> subshapes);

        virtual int countsubshapesrecursively(void);

        // With this the pointer equalities between subshapes in the origin are replicated in this rawshape:
        void replicatelinks(std::shared_ptr<rawshape> origin);

        // Get the mesh info of the shape:
        virtual int getphysicalregion(void);
        virtual std::vector<double>* getcoords(void);
        virtual std::vector<std::vector<int>>* getelems(void);


        virtual std::shared_ptr<rawshape> getpointer(void);

        virtual void mesh(void);
};


#endif
