// sparselizard - Copyright (C) 2017-2018 A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


// This object holds a geometrical shape.

#ifndef RAWEXTRUSION_H
#define RAWEXTRUSION_H

#include <iostream>
#include <vector>

#include "expression.h"

#include "rawshape.h"
#include "rawline.h"
#include "rawpoint.h"
#include "geotools.h"

class rawextrusion: public rawshape
{   

    private:

		int myphysicalregion = -1;
		
		// Son shapes:
		std::vector<std::shared_ptr<rawshape>> sons = {};

		// Coordinates of the nodes in the mesh:
		std::vector<double> mycoords = {};
		// Elements in the mesh:
		std::vector<std::vector<int>> myelems = std::vector<std::vector<int>>(8, std::vector<int>(0));


		// Number of node layers in the extrusion:
		int mynumlayers;

		// Extrusion length:
		double myheight;

		// Dimension of the extruded shape:
		int mydimension;


		// Unextruded shapes:
		std::vector<std::shared_ptr<rawshape>> myunextrudedregions = {};

		// Lines that are the contour of the unextruded shape: 
		std::vector<std::shared_ptr<rawshape>> mycontourregions = {};
        
	public:

		rawextrusion(void) {};

		rawextrusion(int physreg, std::vector<std::shared_ptr<rawshape>> contour, std::vector<std::shared_ptr<rawshape>> innerregions, double height, int numlayers);

		std::shared_ptr<rawshape> duplicate(void);

		void setphysicalregion(int physreg);// + WRITE GETCONTOUR + GETBASE + GETTOP
	
		int getdimension(void);

		std::string getname(void);

		std::vector<std::shared_ptr<rawshape>> getsons(void);

		// Get ALL subshapes (sons are included):
		std::vector<std::shared_ptr<rawshape>> getsubshapes(void);
		
		// Get the mesh info of the shape:
		int getphysicalregion(void);
		std::vector<double>* getcoords(void); 
		std::vector<std::vector<int>>* getelems(void);


		std::shared_ptr<rawshape> getpointer(void);

		void mesh(void);

};

#endif
