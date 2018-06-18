// sparselizard - Copyright (C) 2017-2018 A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


// This object holds a geometrical shape.

#ifndef rawarc_H
#define rawarc_H

#include <iostream>
#include <vector>
#include <memory>

#include "geotools.h"
#include "rawshape.h"
#include "rawpoint.h"
#include "rawextrusion.h"

class rawarc: public rawshape
{   

    private:
        
		int myphysicalregion = -1;
		
		int mynummeshpoints;

		// Son shapes:
		std::vector<std::shared_ptr<rawshape>> sons = {};

		// Coordinates of the nodes in the mesh:
		std::vector<double> mycoords = {};
		// Elements in the mesh:
		std::vector<std::vector<int>> myelems = std::vector<std::vector<int>>(8, std::vector<int>(0));

	public:

		rawarc(void) {};
	
		// Provide to this constructor the point rawshapes at the two ends of the arc and at its center.
		// Order is start point, end point, center point.
		rawarc(int physreg, std::vector<std::shared_ptr<rawshape>> inputpoints, int nummeshpoints);

		std::shared_ptr<rawshape> extrude(int physreg, double height, int numlayers);

		std::shared_ptr<rawshape> duplicate(void);

		// Flip the direction:
		void flip(void);

		void setphysicalregion(int physreg);
	
		int getdimension(void);

		std::string getname(void);

		std::vector<std::shared_ptr<rawshape>> getsons(void);
		
		// Get the mesh info of the shape:
		int getphysicalregion(void);
		std::vector<double>* getcoords(void); 
		std::vector<std::vector<int>>* getelems(void);


		std::shared_ptr<rawshape> getpointer(void);

		void mesh(void); // REWRITE TO WORK IN 3D!!!!!!!!!!!!!!!!!!

};

#endif
