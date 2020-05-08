#include "mesh.h"


mesh::mesh(void)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
}

mesh::mesh(std::string filename, int verbosity, bool legacyreader)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(filename, verbosity, legacyreader);
}

mesh::mesh(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(mergeduplicates, meshfiles, verbosity);
}

mesh::mesh(std::vector<shape> inputshapes, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(inputshapes, verbosity);
}

void mesh::load(std::string name, int verbosity, bool legacyreader)
{
    rawmeshptr->load(name, verbosity, legacyreader);
}

void mesh::load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    rawmeshptr->load(mergeduplicates, meshfiles, verbosity);
}

void mesh::load(std::vector<shape> inputshapes, int verbosity)
{
    rawmeshptr->load(inputshapes, verbosity);
}

void mesh::write(std::string name, int verbosity)
{
    rawmeshptr->write(name, verbosity);
}

void mesh::split(int n)
{
    rawmeshptr->split(n);
}

void mesh::shift(int physreg, double x, double y, double z)
{
    rawmeshptr->shift(physreg, x, y, z);
}

void mesh::shift(double x, double y, double z)
{
    rawmeshptr->shift(x, y, z);
}

void mesh::rotate(int physreg, double ax, double ay, double az)
{
    rawmeshptr->rotate(physreg, ax, ay, az);
}

void mesh::rotate(double ax, double ay, double az)
{
    rawmeshptr->rotate(ax, ay, az);
}

void mesh::scale(int physreg, double x, double y, double z)
{
    rawmeshptr->scale(physreg, x, y, z);
}

void mesh::scale(double x, double y, double z)
{
    rawmeshptr->scale(x, y, z);
}

int mesh::getmeshdimension(void)
{
    return rawmeshptr->getmeshdimension();
}

std::vector<int> mesh::getphysicalregionnumbers(int dim)
{
    return rawmeshptr->getphysicalregionnumbers(dim);
}

void mesh::regionskin(int newphysreg, int physregtoskin)
{
    rawmeshptr->regionskin(newphysreg, physregtoskin);
}

void mesh::boxselection(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit)
{
    rawmeshptr->boxselection(newphysreg, physregtobox, selecteddim, boxlimit);
}

void mesh::sphereselection(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius)
{
    rawmeshptr->sphereselection(newphysreg, physregtosphere, selecteddim, centercoords, radius);
}

void mesh::use(void)
{
    universe::mymesh = rawmeshptr;
}

