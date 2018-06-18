#include "rawdisk.h"


rawdisk::rawdisk(int physreg, std::shared_ptr<rawshape> centerpoint, double radius, int nummeshpts)
{
	if (centerpoint->getdimension() != 0)
	{
		std::cout << "Error in 'rawdisk' object: expected a point shape for the disk center" << std::endl;
		abort();
	}

	if (nummeshpts%4 != 0)
	{
		std::cout << "Error in 'rawdisk' object: the structured disk meshing requires a number of mesh nodes multiple of four" << std::endl;
		abort();
	}

	if (radius == 0)
	{
		std::cout << "Error in 'rawdisk' object: disk radius cannot be zero" << std::endl;
		abort();
	}

	myphysicalregion = physreg;

	mynummeshpoints = nummeshpts;

	mycenterpoint = centerpoint;

	myradius = radius;

	// Sons will be defined when meshing

	mesh();
}


std::shared_ptr<rawshape> rawdisk::extrude(int physreg, double height, int numlayers)
{
	return std::shared_ptr<rawextrusion>(new rawextrusion(physreg, sons, {shared_from_this()}, height, numlayers));
}

std::shared_ptr<rawshape> rawdisk::duplicate(void)
{
    std::shared_ptr<rawdisk> out(new rawdisk);
    *out = *this;
	out->myphysicalregion = -1;

	out->sons = geotools::duplicate(sons);
	out->mycenterpoint = mycenterpoint->duplicate();

    return out;	
}

void rawdisk::setphysicalregion(int physreg)
{
	myphysicalregion = physreg;
}

int rawdisk::getdimension(void)
{
	return 2;
}

std::string rawdisk::getname(void)
{
	return "disk";
}


std::vector<std::shared_ptr<rawshape>> rawdisk::getsons(void) 
{ 
	return sons; 
}

int rawdisk::getphysicalregion(void) 
{ 
	return myphysicalregion; 
}

std::vector<double>* rawdisk::getcoords(void) 
{ 
	return &mycoords; 
}

std::vector<std::vector<int>>* rawdisk::getelems(void) 
{ 
	return &myelems; 
}

std::shared_ptr<rawshape> rawdisk::getpointer(void) 
{ 
	return shared_from_this(); 
}


void rawdisk::mesh(void)
{
	// Get the coordinates of the center node:
	std::vector<double> centercoords = *(mycenterpoint->getcoords());

	if (centercoords[2] != 0.0)
	{
		std::cout << "Error in 'rawdisk' object: disk must be defined in xy plane at zero z" << std::endl;
		abort();
	}
	
	// Define the points for the 4 outer arcs:
	std::shared_ptr<rawshape> p1 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]+myradius, centercoords[1], centercoords[2]}));
	std::shared_ptr<rawshape> p2 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]+myradius, centercoords[2]}));
	std::shared_ptr<rawshape> p3 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]-myradius, centercoords[1], centercoords[2]}));
	std::shared_ptr<rawshape> p4 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]-myradius, centercoords[2]}));

	// Define the points for the 4 corners of the inner square:
	double squareratio = 0.5;
	
	std::shared_ptr<rawshape> p5 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]+squareratio*myradius, centercoords[1], centercoords[2]}));
	std::shared_ptr<rawshape> p6 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]+squareratio*myradius, centercoords[2]}));
	std::shared_ptr<rawshape> p7 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0]-squareratio*myradius, centercoords[1], centercoords[2]}));
	std::shared_ptr<rawshape> p8 = std::shared_ptr<rawpoint>(new rawpoint(-1, {centercoords[0], centercoords[1]-squareratio*myradius, centercoords[2]}));


	// Create the disk contour (4 arcs):
	std::shared_ptr<rawshape> arc1 = std::shared_ptr<rawarc>(new rawarc(-1, {p1,p2}, mynummeshpoints/4));
	std::shared_ptr<rawshape> arc2 = std::shared_ptr<rawarc>(new rawarc(-1, {p2,p3}, mynummeshpoints/4));
	std::shared_ptr<rawshape> arc3 = std::shared_ptr<rawarc>(new rawarc(-1, {p3,p4}, mynummeshpoints/4));
	std::shared_ptr<rawshape> arc4 = std::shared_ptr<rawarc>(new rawarc(-1, {p4,p1}, mynummeshpoints/4));

	// Create the remaining lines:
	std::shared_ptr<rawshape> l1 = std::shared_ptr<rawline>(new rawline(-1, {p1,p5}, mynummeshpoints/4));
	std::shared_ptr<rawshape> l2 = std::shared_ptr<rawline>(new rawline(-1, {p2,p6}, mynummeshpoints/4));
	std::shared_ptr<rawshape> l3 = std::shared_ptr<rawline>(new rawline(-1, {p3,p7}, mynummeshpoints/4));
	std::shared_ptr<rawshape> l4 = std::shared_ptr<rawline>(new rawline(-1, {p4,p8}, mynummeshpoints/4));

	std::shared_ptr<rawshape> l5 = std::shared_ptr<rawline>(new rawline(-1, {p5,p6}, mynummeshpoints/4));
	std::shared_ptr<rawshape> l6 = std::shared_ptr<rawline>(new rawline(-1, {p6,p7}, mynummeshpoints/4));
	std::shared_ptr<rawshape> l7 = std::shared_ptr<rawline>(new rawline(-1, {p7,p8}, mynummeshpoints/4));
	std::shared_ptr<rawshape> l8 = std::shared_ptr<rawline>(new rawline(-1, {p8,p5}, mynummeshpoints/4));

	// Create the 5 quadrangles that make the disk:
	std::shared_ptr<rawshape> q1 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, geotools::orient({arc1,l2,l5,l1})));
	std::shared_ptr<rawshape> q2 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, geotools::orient({arc2,l3,l6,l2})));
	std::shared_ptr<rawshape> q3 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, geotools::orient({arc3,l4,l7,l3})));
	std::shared_ptr<rawshape> q4 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, geotools::orient({arc4,l1,l8,l4})));
	std::shared_ptr<rawshape> q5 = std::shared_ptr<rawquadrangle>(new rawquadrangle(myphysicalregion, geotools::orient({l5,l6,l7,l8})));



	///// Group the mesh of the 5 quadrangles:

	sons = {arc1, arc2, arc3, arc4};

	mycoords = geotools::appendcoords({q1,q2,q3,q4,q5});
	myelems = geotools::appendelems({q1,q2,q3,q4,q5});
}



