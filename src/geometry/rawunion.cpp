#include "rawunion.h"


rawunion::rawunion(int physreg, std::vector<std::shared_ptr<rawshape>> input)
{
    mybuildingblocks = input;
    
    if (input.size() == 0)
    {
        std::cout << "Error in 'rawunion' object: expected at least one shape to unite" << std::endl;
        abort();
    }
    
    int mydimension = input[0]->getdimension();
    for (int i = 1; i < input.size(); i++)
    {
        if (input[i]->getdimension() != mydimension)
        {
            std::cout << "Error in 'rawunion' object: all shapes must have the same dimension" << std::endl;
            abort();
        }
    }

    myphysicalregion = physreg;
        
    mesh();
}

std::shared_ptr<rawshape> rawunion::extrude(int physreg, double height, int numlayers)
{
    std::vector<std::shared_ptr<rawshape>> extrudedshapes(mybuildingblocks.size());
    for (int i = 0; i < mybuildingblocks.size(); i++)
        extrudedshapes[i] = std::shared_ptr<rawextrusion>(new rawextrusion(physreg, mybuildingblocks[i], height, numlayers));

    return std::shared_ptr<rawunion>(new rawunion(physreg, extrudedshapes));
}

std::shared_ptr<rawshape> rawunion::duplicate(void)
{
    std::shared_ptr<rawunion> out(new rawunion);
    *out = *this;

    out->sons = geotools::duplicate(sons);
    out->mybuildingblocks = geotools::duplicate(mybuildingblocks);

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawunion::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawunion::getdimension(void)
{
    return mybuildingblocks[0]->getdimension();
}

std::string rawunion::getname(void)
{
    return "union";
}


std::vector<std::shared_ptr<rawshape>> rawunion::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawunion::getsubshapes(void)
{
    return geotools::concatenate({sons,mybuildingblocks});
}

void rawunion::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    for (int i = 0; i < sons.size(); i++)
        sons[i] = subshapes[i];
        
    for (int i = 0; i < mybuildingblocks.size(); i++)
        mybuildingblocks[i] = subshapes[sons.size()+i];
}

int rawunion::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawunion::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawunion::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawunion::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawunion::mesh(void)
{
    // Append the coordinate of the building blocks:
    mycoords = geotools::appendcoords(mybuildingblocks);
    // Append the elements of the building blocks:
    myelems = geotools::appendelems(mybuildingblocks);
    
    // Create the list of sons:
    std::vector< std::vector<std::shared_ptr<rawshape>> > sonsblocks(mybuildingblocks.size());
    for (int i = 0; i < mybuildingblocks.size(); i++)
        sonsblocks[i] = mybuildingblocks[i]->getsons();
    
    sons = geotools::concatenate(sonsblocks);
}



