#include "rawtriangle.h"


rawtriangle::rawtriangle(int physreg, std::vector<std::shared_ptr<rawshape>> inputpoints, std::vector<int> nummeshpoints)
{
	if (inputpoints.size() != 3 || inputpoints[0]->getdimension() != 0 || inputpoints[1]->getdimension() != 0 || inputpoints[2]->getdimension() != 0)
	{
		std::cout << "Error in 'rawtriangle' object: expected three points in the triangle definition" << std::endl;
		abort();
	}

	if (nummeshpoints.size() != 3)
	{
		std::cout << "Error in 'rawtriangle' object: expected a vector of length three to define the number of mesh nodes" << std::endl;
		abort();
	}

	std::vector<std::shared_ptr<rawshape>> lns(3);
	
	for (int i = 0; i < 3; i++)
		lns[i] = std::shared_ptr<rawline>(new rawline(-1, {inputpoints[i], inputpoints[(i+1)%3]}, nummeshpoints[i]));


	myphysicalregion = physreg;

	sons = lns;

	mesh();
}

rawtriangle::rawtriangle(int physreg, std::vector<std::shared_ptr<rawshape>> inputlines)
{
	if (inputlines.size() != 3 || inputlines[0]->getdimension() != 1 || inputlines[1]->getdimension() != 1 || inputlines[2]->getdimension() != 1)
	{
		std::cout << "Error in 'rawtriangle' object: expected three lines in the triangle definition" << std::endl;
		abort();
	}

	myphysicalregion = physreg;

	sons = geotools::orient(inputlines);

	mesh();
}



std::shared_ptr<rawshape> rawtriangle::extrude(int physreg, double height, int numlayers)
{
	return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, shared_from_this(), height, numlayers));
}

std::shared_ptr<rawshape> rawtriangle::duplicate(void)
{
    std::shared_ptr<rawtriangle> out(new rawtriangle);
    *out = *this;

	out->sons = geotools::duplicate(sons);

	out->replicatelinks(shared_from_this());

    return out;	
}

void rawtriangle::setphysicalregion(int physreg)
{
	myphysicalregion = physreg;
}

int rawtriangle::getdimension(void)
{
	return 2;
}

std::string rawtriangle::getname(void)
{
	return "triangle";
}


std::vector<std::shared_ptr<rawshape>> rawtriangle::getsons(void) 
{ 
	return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawtriangle::getsubshapes(void)
{
	return sons;
}

void rawtriangle::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
	sons = subshapes;
}

int rawtriangle::getphysicalregion(void) 
{ 
	return myphysicalregion; 
}

std::vector<double>* rawtriangle::getcoords(void) 
{ 
	return &mycoords; 
}

std::vector<std::vector<int>>* rawtriangle::getelems(void) 
{ 
	return &myelems; 
}

std::shared_ptr<rawshape> rawtriangle::getpointer(void) 
{ 
	return shared_from_this(); 
}


void rawtriangle::mesh(void)
{
	// Get the node coordinates on the lines:
	std::vector<std::vector<double>*> linescoords(3);
	std::vector<int> nummeshpts(3);

	for (int i = 0; i < 3; i++)
	{
		linescoords[i] = sons[i]->getcoords();
		nummeshpts[i] = linescoords[i]->size()/3;
	}

	// Give an error if the edges do not have the same number of mesh nodes:
	if (nummeshpts[0] != nummeshpts[1] || nummeshpts[0] != nummeshpts[2])
	{
		std::cout << "Error in 'rawtriangle' object: the number of nodes on all edges should be equal" << std::endl;
		abort(); 
	}

	int n = nummeshpts[0];

	// Preallocate the coords and elems containers:
	mycoords.resize(3* 0.5*n*(n+1));
	myelems[2].resize(3*(n-1));
	myelems[3].resize(4* 0.5*(n-2)*(n-1));

	// Get the coordinates of the corner nodes:
	std::vector<double> xcorner = {linescoords[0]->at(0), linescoords[1]->at(0), linescoords[2]->at(0)};
	std::vector<double> ycorner = {linescoords[0]->at(1), linescoords[1]->at(1), linescoords[2]->at(1)};
	std::vector<double> zcorner = {linescoords[0]->at(2), linescoords[1]->at(2), linescoords[2]->at(2)};


	// Loop on all layers in the ki direction:
	int currentnode = 0;
	for (int i = 0; i < n; i++)
	{	
		// Coordinates of the first and last nodes in the line linking 
		// the current node in line 0 and its counterpart in line 1:
		double x1 = linescoords[0]->at(3*i+0);
		double xn = linescoords[1]->at(3*(n-1-i)+0);
		double y1 = linescoords[0]->at(3*i+1);
		double yn = linescoords[1]->at(3*(n-1-i)+1);
		double z1 = linescoords[0]->at(3*i+2);
		double zn = linescoords[1]->at(3*(n-1-i)+2);

		for (int j = 0; j < n-i; j++)
		{

			if (i == 0 && j == n-1)
			{
				mycoords[3*currentnode + 0] = xcorner[2];
				mycoords[3*currentnode + 1] = ycorner[2];
				mycoords[3*currentnode + 2] = zcorner[2];
				currentnode++;
				continue;
			}
			if (i == n-1 && j == 0)
			{
				mycoords[3*currentnode + 0] = xcorner[1];
				mycoords[3*currentnode + 1] = ycorner[1];
				mycoords[3*currentnode + 2] = zcorner[1];
				currentnode++;
				continue;
			}
			if (i == 0 && j == 0)
			{
				mycoords[3*currentnode + 0] = xcorner[0];
				mycoords[3*currentnode + 1] = ycorner[0];
				mycoords[3*currentnode + 2] = zcorner[0];
				currentnode++;
				continue;
			}

			double xh1 = linescoords[2]->at(3*(n-1-j)+0);
			double xhn = linescoords[1]->at(3*j+0);
			double yh1 = linescoords[2]->at(3*(n-1-j)+1);
			double yhn = linescoords[1]->at(3*j+1);
			double zh1 = linescoords[2]->at(3*(n-1-j)+2);
			double zhn = linescoords[1]->at(3*j+2);

			double mu = (double)i/(n-1-j);
			double lambda = (double)j/(n-1-i);

			mycoords[3*currentnode + 0] = 0.5*( (1-lambda)*x1 + lambda*xn + (1-mu)*xh1 + mu*xhn );
			mycoords[3*currentnode + 1] = 0.5*( (1-lambda)*y1 + lambda*yn + (1-mu)*yh1 + mu*yhn );
			mycoords[3*currentnode + 2] = 0.5*( (1-lambda)*z1 + lambda*zn + (1-mu)*zh1 + mu*zhn );

			currentnode++;
		}
	}

	// Add the triangle elements:
	currentnode = 0;
	for (int i = 0; i < n-1; i++)
	{
		currentnode += n-i;

		myelems[2][3*i+0] = currentnode-1;
		myelems[2][3*i+1] = currentnode-2;
		myelems[2][3*i+2] = currentnode+n-(i+1)-1;
	}

	// Add the quadrangle elements:
	int currentelem = 0; currentnode = 0;
	for (int i = 0; i < n-2; i++)
	{
		for (int j = 0; j < n-2-i; j++)
		{
			myelems[3][4*currentelem+0] = currentnode+j+1;
			myelems[3][4*currentelem+1] = currentnode+j;
			myelems[3][4*currentelem+2] = currentnode+n-i+j;
			myelems[3][4*currentelem+3] = currentnode+n-i+j+1;

			currentelem++;
		}

		currentnode += n-i;
	}
}



