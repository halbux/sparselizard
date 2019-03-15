#include "iodata.h"


iodata::iodata(int interpolorder, int geointerpolorder, bool isitscalardata, std::vector<double> timevals)
{
	if (interpolorder <= 0 || geointerpolorder <= 0)
	{
        std::cout << "Error in 'iodata' object: cannot have a negative or zero interpolation order" << std::endl;
        abort();
	}

	myinterpolorder = interpolorder;
	mygeointerpolorder = geointerpolorder;
	isscalardata = isitscalardata;
	mytimevals = timevals;

 	mycoords = std::vector<std::vector<std::vector<densematrix>>>(3, std::vector<std::vector<densematrix>>(8, std::vector<densematrix>(0)));
	myscalardata = std::vector<std::vector<densematrix>>(8, std::vector<densematrix>(0));
	myvectordata = std::vector<std::vector<std::vector<densematrix>>>(3, std::vector<std::vector<densematrix>>(8, std::vector<densematrix>(0)));
}

void iodata::combine(void)
{
	// Loop on every element type:
	for (int i = 0; i < 8; i++)
	{
		// Skip if there is a single block:
		if (mycoords[0][i].size() <= 1)
			continue;
	
		for (int s = 0; s < 3; s++)
			mycoords[s][i] = {densematrix(mycoords[s][i])};
		myscalardata[i] = {densematrix(myscalardata[i])};
		for (int comp = 0; comp < 3; comp++)
			myvectordata[comp][i] = {densematrix(myvectordata[comp][i])};
	}
}

bool iodata::isscalar(void) { return isscalardata; }
int iodata::getinterpolorder(void) { return myinterpolorder; };
int iodata::getgeointerpolorder(void) { return mygeointerpolorder; };

std::vector<double> iodata::gettimetags(void) { return mytimevals; };

bool iodata::ispopulated(int elemtypenum)
{
	return (mycoords[0][elemtypenum].size() > 0);
}

void iodata::addcoordinates(int elemtypenum, densematrix xcoords, densematrix ycoords, densematrix zcoords)
{
	mycoords[0][elemtypenum].push_back(xcoords);
	mycoords[1][elemtypenum].push_back(ycoords);
	mycoords[2][elemtypenum].push_back(zcoords);
}

void iodata::addscalardata(int elemtypenum, densematrix vals)
{
	if (isscalardata == false)
	{
        std::cout << "Error in 'iodata' object: expected scalar data" << std::endl;
        abort();
	}

	myscalardata[elemtypenum].push_back(vals);
}

void iodata::addvectordata(int elemtypenum, densematrix compxvals, densematrix compyvals, densematrix compzvals)
{
	if (isscalardata == true)
	{
        std::cout << "Error in 'iodata' object: expected vector data" << std::endl;
        abort();
	}

	myvectordata[0][elemtypenum].push_back(compxvals);
	myvectordata[1][elemtypenum].push_back(compyvals);
	myvectordata[2][elemtypenum].push_back(compzvals);
}

std::vector<densematrix> iodata::getcoordinates(int elemtypenum)
{
	combine();
	return {mycoords[0][elemtypenum][0], mycoords[1][elemtypenum][0], mycoords[2][elemtypenum][0]};
}

densematrix iodata::getscalardata(int elemtypenum)
{
	combine();
	return myscalardata[elemtypenum][0];
}

std::vector<densematrix> iodata::getvectordata(int elemtypenum)
{
	combine();
	return {myvectordata[0][elemtypenum][0], myvectordata[1][elemtypenum][0], myvectordata[2][elemtypenum][0]};
}

