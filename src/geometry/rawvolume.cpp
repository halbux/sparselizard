#include "rawvolume.h"


rawvolume::rawvolume(int physreg, std::vector<double>& allcoords, std::vector<std::vector<int>>& allelems)
{
    myphysicalregion = physreg;

    mycoords = allcoords;
    myelems = allelems;
}



std::shared_ptr<rawshape> rawvolume::duplicate(void)
{
    std::shared_ptr<rawvolume> out(new rawvolume);
    *out = *this;

    out->sons = geotools::duplicate(sons);

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawvolume::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawvolume::getdimension(void)
{
    return 3;
}

std::string rawvolume::getname(void)
{
    return "volume";
}


std::vector<std::shared_ptr<rawshape>> rawvolume::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawvolume::getsubshapes(void)
{
    return sons;
}

void rawvolume::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    sons = subshapes;
}

int rawvolume::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawvolume::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawvolume::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawvolume::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawvolume::mesh(void)
{

}



