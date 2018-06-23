#include "geotools.h"


double geotools::acos(double arg)
{
    double pi = 3.1415926535897932384;

	if (arg >= 1)
		return 0.0;
	if (arg <= -1)
		return pi;
	if (std::abs(arg) < 1)
		return std::acos(arg);
}

std::vector<std::shared_ptr<rawshape>> geotools::coordstopoints(std::vector<double> coords)
{
	if (coords.size()%3 == 0)
	{
		int numpts = coords.size()/3;	

		std::vector<std::shared_ptr<rawshape>> pts(numpts);
		for (int i = 0; i < numpts; i++)
			pts[i] = std::shared_ptr<rawpoint>(new rawpoint(-1, {coords[3*i+0], coords[3*i+1], coords[3*i+2]}));

		return pts;
	}
	else
	{
		std::cout << "Error in 'geotool' namespace: length of coordinate vector should be a multiple of three" << std::endl;
		abort();
	}
}

double geotools::getdistance(std::vector<double> pt1coords, std::vector<double> pt2coords)
{
	return std::sqrt( std::pow(pt1coords[0]-pt2coords[0],2) + std::pow(pt1coords[1]-pt2coords[1],2) + std::pow(pt1coords[2]-pt2coords[2],2) );
}

std::vector<double> geotools::getplaneangles(std::vector<double> p1, std::vector<double> p2, std::vector<double> p3)
{
    double pi = 3.1415926535897932384;

	// Roundoff noise threshold:
	double threshold = 1e-10;

	// Vector between two points:
	std::vector<double> v12 = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
	std::vector<double> v13 = {p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]};
	double v12norm = std::sqrt(v12[0]*v12[0]+v12[1]*v12[1]+v12[2]*v12[2]);

	// Shortcut if all points are in a plane parallel to the xy plane.
	// This allows colinear points in 2D.
	if (std::abs(v12[2])/v12norm < threshold && std::abs(v13[2])/v12norm < threshold)
		return std::vector<double>{0.0,0.0};

	// Normal [a b c] to the plane is the cross product:
	std::vector<double> vnormal = {v12[1]*v13[2]-v12[2]*v13[1], v12[2]*v13[0]-v12[0]*v13[2], v12[0]*v13[1]-v12[1]*v13[0]};
	double normalnorm = std::sqrt(vnormal[0]*vnormal[0]+vnormal[1]*vnormal[1]+vnormal[2]*vnormal[2]);

	// If in the xz plane:
	if (std::abs(vnormal[0])/normalnorm < threshold && std::abs(vnormal[2])/normalnorm < threshold)
		return std::vector<double>{0.0,90.0};
	// If in the yz plane:
	if (std::abs(vnormal[1])/normalnorm < threshold && std::abs(vnormal[2])/normalnorm < threshold)
		return std::vector<double>{90.0,0.0};

	// Plane equation is ax + by + cz = d.
	// Setting x then y to zero gives the angles we need:
	return std::vector<double>{360/2/pi*std::atan(-vnormal[0]/vnormal[2]), 360/2/pi*std::atan(-vnormal[1]/vnormal[2])};
}

std::vector<double> geotools::flipcoords(std::vector<double>& input)
{
	int numnodes = input.size()/3;

	std::vector<double> output(3*numnodes);

	for (int i = 0; i < numnodes; i++)
	{
		output[3*i+0] = input[3*(numnodes-1-i)+0];
		output[3*i+1] = input[3*(numnodes-1-i)+1];
		output[3*i+2] = input[3*(numnodes-1-i)+2];
	}
	
	return output;
}

std::vector<std::shared_ptr<rawshape>> geotools::orient(std::vector<std::shared_ptr<rawshape>> input)
{
	if (input.size() < 2)
		return input;

	// Get the coordinates of the points in the two first lines:
	std::vector<double> p1coord = *(input[0]->getsons()[0]->getcoords());
	std::vector<double> p2coord = *(input[0]->getsons()[1]->getcoords());
	std::vector<double> p3coord = *(input[1]->getsons()[0]->getcoords());
	std::vector<double> p4coord = *(input[1]->getsons()[1]->getcoords());

	// Find which pair between the two lines is the closest in distance:
	double d13 = std::sqrt( std::pow(p1coord[0]-p3coord[0],2) + std::pow(p1coord[1]-p3coord[1],2) + std::pow(p1coord[2]-p3coord[2],2) );
	double d14 = std::sqrt( std::pow(p1coord[0]-p4coord[0],2) + std::pow(p1coord[1]-p4coord[1],2) + std::pow(p1coord[2]-p4coord[2],2) );
	double d23 = std::sqrt( std::pow(p2coord[0]-p3coord[0],2) + std::pow(p2coord[1]-p3coord[1],2) + std::pow(p2coord[2]-p3coord[2],2) );
	double d24 = std::sqrt( std::pow(p2coord[0]-p4coord[0],2) + std::pow(p2coord[1]-p4coord[1],2) + std::pow(p2coord[2]-p4coord[2],2) );

	// This will tell which lines must be flipped:
	std::vector<bool> flipline(input.size(), false);

	if (d13 < d14 && d13 < d23 && d13 < d24)
		flipline[0] = true;
	if (d14 < d13 && d14 < d23 && d14 < d24)
		flipline[0] = flipline[1] = true;
	if (d24 < d13 && d24 < d23 && d24 < d14)
		flipline[1] = true;

	
	std::vector<double> prevnodecoord;
	if (flipline[1] == false)
		prevnodecoord = p4coord;
	else
		prevnodecoord = p3coord;

	for (int i = 2; i < input.size(); i++)
	{
		// Get the coordinates of the points in the current line:
		std::vector<double> p1c = *(input[i]->getsons()[0]->getcoords());
		std::vector<double> p2c = *(input[i]->getsons()[1]->getcoords());

		// Find which point is closest to the previous one:
		double d1 = std::sqrt( std::pow(p1c[0]-prevnodecoord[0],2) + std::pow(p1c[1]-prevnodecoord[1],2) + std::pow(p1c[2]-prevnodecoord[2],2) );
		double d2 = std::sqrt( std::pow(p2c[0]-prevnodecoord[0],2) + std::pow(p2c[1]-prevnodecoord[1],2) + std::pow(p2c[2]-prevnodecoord[2],2) );
			
		if (d1 < d2)
			prevnodecoord = p2c;
		else
		{
			flipline[i] = true;
			prevnodecoord = p1c;
		}
	}


	// Flip the lines that must be flipped:
	std::vector<std::shared_ptr<rawshape>> output = geotools::duplicate(input);

	for (int i = 0; i < output.size(); i++)
	{
		if (flipline[i] == true)
			output[i]->flip();
	}

	return output;
}

std::vector< std::shared_ptr<rawshape> > geotools::getrawshapes(std::vector<shape> shapes)
{
	std::vector< std::shared_ptr<rawshape> > rawshapes(shapes.size());

	for (int i = 0; i < shapes.size(); i++)
	{
		std::shared_ptr<rawshape> currentrawshapeptr = shapes[i].getpointer();

		if (currentrawshapeptr != NULL)	
			rawshapes[i] = currentrawshapeptr;
		else
		{
			std::cout << "Error in 'geotool' namespace: encountered an undefined shape (NULL rawshape pointer)" << std::endl;
			abort();
		}
	}
	
	return rawshapes;
}

std::vector<shape> geotools::getshapes(std::vector< std::shared_ptr<rawshape> > rawshapes)
{
	std::vector<shape> shapes(rawshapes.size());

	for (int i = 0; i < rawshapes.size(); i++)
		shapes[i] = shape(rawshapes[i]);
	
	return shapes;
}

std::vector< std::shared_ptr<rawshape> > geotools::flip(std::vector< std::shared_ptr<rawshape> > input)
{
	std::vector< std::shared_ptr<rawshape> > output(input.size());

	for (int i = 0; i < input.size(); i++)
		output[i] = input[input.size()-1-i];	

	return output;
}

std::vector<std::shared_ptr<rawshape>> geotools::duplicate(std::vector<std::shared_ptr<rawshape>> input)
{
	std::vector<std::shared_ptr<rawshape>> output(input.size());

	for (int i = 0; i < input.size(); i++)
		output[i] = input[i]->duplicate();

	return output;
}

std::vector<std::shared_ptr<rawshape>> geotools::concatenate(std::vector<std::vector<std::shared_ptr<rawshape>>> input)
{
	// Compute the total number of rawshapes:
	int numrawshapes = 0;
	for (int i = 0; i < input.size(); i++)
		numrawshapes += input[i].size();

	std::vector<std::shared_ptr<rawshape>> output(numrawshapes);

	int index = 0;
	for (int i = 0; i < input.size(); i++)
	{
		for (int j = 0; j < input[i].size(); j++)
		{
			output[index] = input[i][j];
			index++;
		}
	}

	return output;
}

std::vector<double> geotools::appendcoords(std::vector<std::shared_ptr<rawshape>> rawshapes)
{
	std::vector< std::vector<double>* > coordsptrs(rawshapes.size());

	for (int i = 0; i < rawshapes.size(); i++)
		coordsptrs[i] = rawshapes[i]->getcoords();
	

	// First calculate the total size:
	int totallen = 0;
	for (int i = 0; i < coordsptrs.size(); i++)
		totallen += coordsptrs[i]->size();

	std::vector<double> output(totallen);

	int index = 0;
	for (int i = 0; i < coordsptrs.size(); i++)
	{
		for (int j = 0; j < coordsptrs[i]->size(); j++)
		{
			output[index] = coordsptrs[i]->at(j);
			index++;
		}
	}
	return output;
}

std::vector<std::vector<int>> geotools::appendelems(std::vector<std::shared_ptr<rawshape>> rawshapes)
{
	std::vector< std::vector<std::vector<int>>* > elemsptrs(rawshapes.size());

	for (int i = 0; i < rawshapes.size(); i++)
		elemsptrs[i] = rawshapes[i]->getelems();
	

	std::vector<std::vector<int>> output(8);

	// Loop on all elements:
	for (int e = 0; e < 8; e++)
	{
		// First calculate the total size:
		int totallen = 0;
		for (int i = 0; i < elemsptrs.size(); i++)
			totallen += elemsptrs[i]->at(e).size();

		output[e] = std::vector<int>(totallen);

		int index = 0;
		for (int i = 0; i < elemsptrs.size(); i++)
		{
			for (int j = 0; j < elemsptrs[i]->at(e).size(); j++)
			{
				output[e][index] = elemsptrs[i]->at(e)[j];
				index++;
			}
		}
	}
	return output;
}



