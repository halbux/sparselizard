// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This namespace provides tools for the geometry operations.

#ifndef GEOTOOLS_H
#define GEOTOOLS_H

#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <memory>
#include "shape.h"
#include "rawshape.h"
#include "rawpoint.h"

namespace geotools
{
    // An implementation of acos that is safe with respect to roundoff noise on the argument:
    double acos(double arg);

    // Convert a list of point coordinates into point objects:
    std::vector<std::shared_ptr<rawshape>> coordstopoints(std::vector<double> coords);

    // Get the distance between two points:
    double getdistance(std::vector<double> pt1coords, std::vector<double> pt2coords);
    double getdistance(int pt1, int pt2, std::vector<double>& coords);

    // Get the angle (in degrees) which rotates the plane defined by the 3 input points
    // around the x (or y) axis to bring it parallel to the y axis (or x axis).
    // Get the angle around the x axis with "xrot" and around the y axis with "yrot"
    double getplanerotation(std::string xy, std::vector<double> p1, std::vector<double> p2, std::vector<double> p3);

    // Rotate a vector of coordinates by alphax, alphay and alphaz degrees around the x, y and z axis respectively:
    void rotate(double alphax, double alphay, double alphaz, std::vector<double>* coords);

    // Flip the nodes in a coordinate vector (3 coordinates per node):
    std::vector<double> flipcoords(std::vector<double>& input);

    // Orient a list of line shapes to have them all pointing to the next one:
    std::vector<std::shared_ptr<rawshape>> orient(std::vector<std::shared_ptr<rawshape>> input);

    // Transform a vector of shapes into a vector of rawshapes:
    std::vector< std::shared_ptr<rawshape> > getrawshapes(std::vector<shape> shapes);

    // Transform a vector of rawshapes into a vector of shapes:
    std::vector<shape> getshapes(std::vector< std::shared_ptr<rawshape> > rawshapes);

    // Transform a vector of rawshape shared pointers to a vector of rawshape pointers:
    std::vector<rawshape*> getpointers(std::vector< std::shared_ptr<rawshape> > sharedptrs);

    // Flip the rawshape vector direction:
    std::vector< std::shared_ptr<rawshape> > flip(std::vector< std::shared_ptr<rawshape> > input);

    // Unique a list of rawshape pointers:
    std::vector<rawshape*> unique(std::vector<rawshape*> ptrs);

    // Duplicate a list of rawshapes:
    std::vector<std::shared_ptr<rawshape>> duplicate(std::vector<std::shared_ptr<rawshape>> input);

    // Concatenate lists of rawshapes:
    std::vector<std::shared_ptr<rawshape>> concatenate(std::vector<std::vector<std::shared_ptr<rawshape>>> input);

    // Get an int vector to sort a vector of rawshape pointers:
    void sortrawshapepointers(std::vector<rawshape*>& tosort, std::vector<int>& reorderingvector);

    // Append the coordinate of multiple rawshapes:
    std::vector<double> appendcoords(std::vector<std::shared_ptr<rawshape>> rawshapes);
    // Append the elements of multiple rawshapes:
    std::vector<std::vector<int>> appendelems(std::vector<std::shared_ptr<rawshape>> rawshapes);
};

#endif
