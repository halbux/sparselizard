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
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

mesh::mesh(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(mergeduplicates, meshfiles, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

mesh::mesh(std::vector<shape> inputshapes, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(inputshapes, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

void mesh::load(std::string name, int verbosity, bool legacyreader)
{
    errorifloaded();
    rawmeshptr->load(name, verbosity, legacyreader);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

void mesh::load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    errorifloaded();
    rawmeshptr->load(mergeduplicates, meshfiles, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

void mesh::load(std::vector<shape> inputshapes, int verbosity)
{
    errorifloaded();
    rawmeshptr->load(inputshapes, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

void mesh::write(std::string name, int verbosity)
{
    rawmeshptr->gethadaptedpointer()->write(name, verbosity);
}

void mesh::setadaptivity(expression criterion, int lownumsplits, int highnumsplits)
{
    if (not(criterion.isscalar()))
    {
        std::cout << "Error in 'mesh' object: expected a scalar criterion for h-adaptivity" << std::endl;
        abort();   
    }
    // The criterion cannot be multiharmonic:
    std::vector<int> alldisjregs(universe::mymesh->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(criterion.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'mesh' object: cannot have a multiharmonic criterion for h-adaptivity" << std::endl;
        abort();
    }
    
    if (lownumsplits < 0)
    {
        std::cout << "Error in 'mesh' object: in 'setadaptivity' cannot use negative minimum number of splits " << lownumsplits << std::endl;
        abort();   
    }
    if (highnumsplits < lownumsplits)
    {
        std::cout << "Error in 'mesh' object: in 'setadaptivity' the minimum number of splits cannot be larger than the maximum" << std::endl;
        abort();   
    }
    
    rawmeshptr->setadaptivity(criterion, lownumsplits, highnumsplits);
}

void mesh::split(int n)
{
    errorifloaded();
    rawmeshptr->split(n);
}

void mesh::move(int physreg, expression u)
{
    if (rawmeshptr->getmeshnumber() == 0)
        rawmeshptr->move(physreg, u);
    rawmeshptr->gethadaptedpointer()->move(physreg, u);
}

void mesh::move(expression u)
{
    if (rawmeshptr->getmeshnumber() == 0)
        rawmeshptr->move(-1, u);
    rawmeshptr->gethadaptedpointer()->move(-1, u);
}
        
void mesh::shift(int physreg, double x, double y, double z)
{
    rawmeshptr->shift(physreg, x, y, z);
    rawmeshptr->gethadaptedpointer()->shift(physreg, x, y, z);
}

void mesh::shift(double x, double y, double z)
{
    rawmeshptr->shift(-1, x, y, z);
    rawmeshptr->gethadaptedpointer()->shift(-1, x, y, z);
}

void mesh::rotate(int physreg, double ax, double ay, double az)
{
    rawmeshptr->rotate(physreg, ax, ay, az);
    rawmeshptr->gethadaptedpointer()->rotate(physreg, ax, ay, az);
}

void mesh::rotate(double ax, double ay, double az)
{
    rawmeshptr->rotate(-1, ax, ay, az);
    rawmeshptr->gethadaptedpointer()->rotate(-1, ax, ay, az);
}

void mesh::scale(int physreg, double x, double y, double z)
{
    rawmeshptr->scale(physreg, x, y, z);
    rawmeshptr->gethadaptedpointer()->scale(physreg, x, y, z);
}

void mesh::scale(double x, double y, double z)
{
    rawmeshptr->scale(-1, x, y, z);
    rawmeshptr->gethadaptedpointer()->scale(-1, x, y, z);
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

void mesh::layerselection(int newphysreg, int physregtoselectfrom, int physregtostartgrowth, int numlayers)
{
    errorifloaded();
    rawmeshptr->layerselection(newphysreg, physregtoselectfrom, physregtostartgrowth, numlayers);
}

void mesh::regionexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude)
{
    errorifloaded();
    rawmeshptr->regionexclusion(newphysreg, physregtoexcludefrom, physregstoexclude);
}

void mesh::use(void)
{
    universe::mymesh = rawmeshptr->gethadaptedpointer();
}

