#include "mesh.h"


void mesh::errorifloaded(void)
{
    if (isloaded)
    {
        std::cout << "Error in 'mesh' object: cannot perform the requested operation (mesh is already loaded)" << std::endl;
        abort();
    }
}

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
    isloaded = true;
}

mesh::mesh(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(mergeduplicates, meshfiles, verbosity);
    isloaded = true;
}

mesh::mesh(std::vector<shape> inputshapes, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(inputshapes, verbosity);
    isloaded = true;
}

void mesh::load(std::string name, int verbosity, bool legacyreader)
{
    errorifloaded();
    rawmeshptr->load(name, verbosity, legacyreader);
    isloaded = true;
}

void mesh::load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    errorifloaded();
    rawmeshptr->load(mergeduplicates, meshfiles, verbosity);
    isloaded = true;
}

void mesh::load(std::vector<shape> inputshapes, int verbosity)
{
    errorifloaded();
    rawmeshptr->load(inputshapes, verbosity);
    isloaded = true;
}

void mesh::write(std::string name, int verbosity)
{
    if (rawmeshptr->gethadaptedpointer() != NULL)
        rawmeshptr->gethadaptedpointer()->write(name, verbosity);
    else
        rawmeshptr->write(name, verbosity);
}

bool mesh::adapt(int verbosity)
{
    return rawmeshptr->adapth(verbosity);
}

void mesh::setadaptivity(expression criterion, std::vector<field> triggers, int lownumsplits, int highnumsplits, double thresdown, double thresup, double mincritrange)
{
    rawmeshptr->setadaptivity(criterion, triggers, lownumsplits, highnumsplits, thresdown, thresup, mincritrange);
}

void mesh::setadaptivity(expression criterion, std::vector<field> triggers, std::vector<double> thresholds, std::vector<int> numsplits, double thresdown, double thresup, double mincritrange)
{
    rawmeshptr->setadaptivity(criterion, triggers, thresholds, numsplits, thresdown, thresup, mincritrange);
}

void mesh::split(int n)
{
    errorifloaded();
    rawmeshptr->split(n);
}

void mesh::shift(int physreg, double x, double y, double z)
{
    rawmeshptr->shift(physreg, x, y, z);
    rawmeshptr->gethadaptedpointer()->shift(physreg, x, y, z);
}

void mesh::shift(double x, double y, double z)
{
    rawmeshptr->shift(x, y, z);
    rawmeshptr->gethadaptedpointer()->shift(x, y, z);
}

void mesh::rotate(int physreg, double ax, double ay, double az)
{
    rawmeshptr->rotate(physreg, ax, ay, az);
    rawmeshptr->gethadaptedpointer()->rotate(physreg, ax, ay, az);
}

void mesh::rotate(double ax, double ay, double az)
{
    rawmeshptr->rotate(ax, ay, az);
    rawmeshptr->gethadaptedpointer()->rotate(ax, ay, az);
}

void mesh::scale(int physreg, double x, double y, double z)
{
    rawmeshptr->scale(physreg, x, y, z);
    rawmeshptr->gethadaptedpointer()->scale(physreg, x, y, z);
}

void mesh::scale(double x, double y, double z)
{
    rawmeshptr->scale(x, y, z);
    rawmeshptr->gethadaptedpointer()->scale(x, y, z);
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
    errorifloaded();
    rawmeshptr->regionskin(newphysreg, physregtoskin);
}

void mesh::boxselection(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit)
{
    errorifloaded();
    rawmeshptr->boxselection(newphysreg, physregtobox, selecteddim, boxlimit);
}

void mesh::sphereselection(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius)
{
    errorifloaded();
    rawmeshptr->sphereselection(newphysreg, physregtosphere, selecteddim, centercoords, radius);
}

void mesh::regionexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude)
{
    errorifloaded();
    rawmeshptr->regionexclusion(newphysreg, physregtoexcludefrom, physregstoexclude);
}

void mesh::use(void)
{
    if (rawmeshptr->gethadaptedpointer() != NULL)
        universe::mymesh = rawmeshptr->gethadaptedpointer();
    else
        universe::mymesh = rawmeshptr;
}

