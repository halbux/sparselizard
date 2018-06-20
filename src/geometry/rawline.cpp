#include "rawline.h"


rawline::rawline(int physreg, std::vector<double> allcoords)
{
	myphysicalregion = physreg;

	mynummeshpoints = allcoords.size()/3;
	int numelems = mynummeshpoints-1;

	if (mynummeshpoints < 2)
	{
		std::cout << "Error in 'rawline' object: expected at least two mesh nodes in the line" << std::endl;
		abort();
	}

	// First and last coordinates correspond to the two corner points:
	std::shared_ptr<rawpoint> p1(new rawpoint(-1, {allcoords[0], allcoords[1], allcoords[2]}));
	std::shared_ptr<rawpoint> p2(new rawpoint(-1, {allcoords[3*(mynummeshpoints-1)+0], allcoords[3*(mynummeshpoints-1)+1], allcoords[3*(mynummeshpoints-1)+2]}));

	sons = {p1, p2};

	mycoords = allcoords;

	// Here the mesh is provided by the user:
	myelems[1].resize(2*numelems);
	for (int i = 0; i < numelems; i++)
	{
		myelems[1][2*i+0] = i;
		myelems[1][2*i+1] = i+1;
	}
}

rawline::rawline(int physreg, std::vector<std::shared_ptr<rawshape>> inputpoints, int nummeshpoints)
{
	if (inputpoints.size() != 2 || inputpoints[0]->getdimension() != 0 || inputpoints[1]->getdimension() != 0)
	{
		std::cout << "Error in 'rawline' object: expected two points in the line definition" << std::endl;
		abort();
	}

	if (nummeshpoints < 2)
	{
		std::cout << "Error in 'rawline' object: expected at least two mesh nodes in the line" << std::endl;
		abort();
	}

	myphysicalregion = physreg;

	mynummeshpoints = nummeshpoints;

	sons = inputpoints;

	mesh();
}
 

std::shared_ptr<rawshape> rawline::extrude(int physreg, double height, int numlayers)
{
	return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, sons, {shared_from_this()}, height, numlayers));
}

std::shared_ptr<rawshape> rawline::duplicate(void)
{
    std::shared_ptr<rawline> out(new rawline);
    *out = *this;

	out->sons = geotools::duplicate(sons);

    return out;	
}

void rawline::flip(void)
{
	sons = {sons[1],sons[0]};

	int numlines = myelems[1].size()/2;

	std::vector<int> flippedlines(2*numlines);

	for (int i = 0; i < 2*numlines; i++)
		flippedlines[i] = myelems[1][2*numlines-1-i];
		
	myelems[1] = flippedlines;
}

void rawline::setphysicalregion(int physreg)
{
	myphysicalregion = physreg;
}

int rawline::getdimension(void)
{
	return 1;
}

std::string rawline::getname(void)
{
	return "line";
}


std::vector<std::shared_ptr<rawshape>> rawline::getsons(void) 
{ 
	return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawline::getsubshapes(void)
{
	return sons;
}

int rawline::getphysicalregion(void) 
{ 
	return myphysicalregion; 
}

std::vector<double>* rawline::getcoords(void) 
{ 
	return &mycoords; 
}

std::vector<std::vector<int>>* rawline::getelems(void) 
{ 
	return &myelems; 
}

std::shared_ptr<rawshape> rawline::getpointer(void) 
{ 
	return shared_from_this(); 
}


void rawline::mesh(void)
{
	int numelems = mynummeshpoints-1;

	mycoords.resize(3*mynummeshpoints);
	myelems[1].resize(2*numelems);

	std::vector<double>* p1coords = sons[0]->getcoords();
	std::vector<double>* p2coords = sons[1]->getcoords();

	for (int i = 0; i < mynummeshpoints; i++)
	{
		// Weight:
		double w = (double)i/(mynummeshpoints-1);
	
		mycoords[3*i+0] = (1-w)*p1coords->at(0) + w*p2coords->at(0);
		mycoords[3*i+1] = (1-w)*p1coords->at(1) + w*p2coords->at(1);
		mycoords[3*i+2] = (1-w)*p1coords->at(2) + w*p2coords->at(2);
	}

	for (int i = 0; i < numelems; i++)
	{
		myelems[1][2*i+0] = i;
		myelems[1][2*i+1] = i+1;
	}
}

