// sparselizard - Copyright (C) 2017-2018 A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


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

		virtual void deform(expression xdeform, expression ydeform, expression zdeform);

		virtual void shift(double shiftx, double shifty, double shiftz);
		virtual void rotate(double alphax, double alphay, double alphaz);


		virtual std::shared_ptr<rawshape> extrude(int physreg, double height, int numlayers);

		virtual std::shared_ptr<rawshape> duplicate(void);

		// Flip a line direction:
		virtual void flip(void);

		virtual void setphysicalregion(int physreg);

		virtual int getdimension(void);

		virtual std::string getname(void);

		virtual std::vector<std::shared_ptr<rawshape>> getsons(void);

		// Get ALL subshapes (sons are included):
		virtual std::vector<std::shared_ptr<rawshape>> getsubshapes(void);

		// Get the mesh info of the shape:
		virtual int getphysicalregion(void);
		virtual std::vector<double>* getcoords(void);
		virtual std::vector<std::vector<int>>* getelems(void);


		virtual std::shared_ptr<rawshape> getpointer(void);

		virtual void mesh(void);
};


#endif
