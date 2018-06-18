// sparselizard - Copyright (C) 2017-2018 A. Halbach
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


// This namespace provides tools for the geometry operations.

#ifndef GEOTOOLS_H
#define GEOTOOLS_H

#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <memory>
#include "shape.h"
#include "rawshape.h"
#include "rawpoint.h"

namespace geotools
{
	// Convert a list of point coordinates into point objects:
	std::vector<std::shared_ptr<rawshape>> coordstopoints(std::vector<double> coords);

	// Orient a list of line shapes to have them all pointing to the next one:
	std::vector<std::shared_ptr<rawshape>> orient(std::vector<std::shared_ptr<rawshape>> input);

	// Transform a vector of shapes into a vector of rawshapes:
	std::vector< std::shared_ptr<rawshape> > getrawshapes(std::vector<shape> shapes);

	// Transform a vector of rawshapes into a vector of shapes:
	std::vector<shape> getshapes(std::vector< std::shared_ptr<rawshape> > rawshapes);

	// Flip the rawshape vector direction:
	std::vector< std::shared_ptr<rawshape> > flip(std::vector< std::shared_ptr<rawshape> > input);

	// Duplicate a list of rawshapes:
	std::vector<std::shared_ptr<rawshape>> duplicate(std::vector<std::shared_ptr<rawshape>> input);

	// Concatenate lists of rawshapes:
	std::vector<std::shared_ptr<rawshape>> concatenate(std::vector<std::vector<std::shared_ptr<rawshape>>> input);

	// Append the coordinate of multiple rawshapes:
    std::vector<double> appendcoords(std::vector<std::shared_ptr<rawshape>> rawshapes);
	// Append the elements of multiple rawshapes:
	std::vector<std::vector<int>> appendelems(std::vector<std::shared_ptr<rawshape>> rawshapes);
};

#endif
