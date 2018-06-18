#include "geotools.h"


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



