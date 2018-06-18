// sparselizard - Copyright (C) 2017-2018 A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


// This object holds a geometrical shape.

#ifndef RAWPOINT_H
#define RAWPOINT_H

#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>

#include "rawshape.h"
#include "rawextrusion.h"

class rawpoint: public rawshape
{   

    private:
        
		int myphysicalregion = -1;

		// Son shapes:
		std::vector<std::shared_ptr<rawshape>> sons = {};

		// Coordinates of the nodes in the mesh:
		std::vector<double> mycoords = {};
		// Elements in the mesh:
		std::vector<std::vector<int>> myelems = std::vector<std::vector<int>>(8, std::vector<int>(0));

	public:
 
		rawpoint(void) {};

		rawpoint(int physreg, std::vector<double> coords);
      
		std::shared_ptr<rawshape> extrude(int physreg, double height, int numlayers);

		std::shared_ptr<rawshape> duplicate(void);

		void setphysicalregion(int physreg);
	
		int getdimension(void);

		std::string getname(void);

		std::vector<std::shared_ptr<rawshape>> getsons(void);
		
		// Get the mesh info of the shape:
		int getphysicalregion(void);
		std::vector<double>* getcoords(void); 
		std::vector<std::vector<int>>* getelems(void);


		std::shared_ptr<rawshape> getpointer(void);
};

#endif
