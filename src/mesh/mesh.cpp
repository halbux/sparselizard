#include "mesh.h"


void mesh::errorifloaded(void)
{
    if (isloaded)
    {
        std::cout << "Error in 'mesh' object: cannot perform the requested operation (mesh is already loaded)" << std::endl;
        abort();
    }
}

void mesh::errorifnotloaded(void)
{
    if (not(isloaded))
    {
        std::cout << "Error in 'mesh' object: cannot perform the requested operation (mesh is not loaded)" << std::endl;
        abort();
    }
}

mesh::mesh(void)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
}

mesh::mesh(std::string filename, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(filename, -1, -1, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

mesh::mesh(std::string filename, int globalgeometryskin, int numoverlaplayers, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(filename, globalgeometryskin, numoverlaplayers, verbosity);
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
    rawmeshptr->load(inputshapes, -1, -1, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

mesh::mesh(std::vector<shape> inputshapes, int globalgeometryskin, int numoverlaplayers, int verbosity)
{
    rawmeshptr = std::shared_ptr<rawmesh>(new rawmesh());
    universe::mymesh = rawmeshptr;
    rawmeshptr->load(inputshapes, globalgeometryskin, numoverlaplayers, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

void mesh::load(std::string name, int verbosity)
{
    load(name, -1, -1, verbosity);
}

void mesh::load(std::string name, int globalgeometryskin, int numoverlaplayers, int verbosity)
{
    errorifloaded();
    rawmeshptr->load(name, globalgeometryskin, numoverlaplayers, verbosity);
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
    load(inputshapes, -1, -1, verbosity);
}

void mesh::load(std::vector<shape> inputshapes, int globalgeometryskin, int numoverlaplayers, int verbosity)
{
    errorifloaded();
    rawmeshptr->load(inputshapes, globalgeometryskin, numoverlaplayers, verbosity);
    universe::mymesh = rawmeshptr->gethadaptedpointer();
    isloaded = true;
}

void mesh::write(std::string name, int verbosity)
{
    errorifnotloaded();
    rawmeshptr->gethadaptedpointer()->write(name, verbosity);
}

void mesh::setadaptivity(expression criterion, int lownumsplits, int highnumsplits)
{
    errorifnotloaded();
    
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
    
    rawmeshptr->setadaptivity(criterion, lownumsplits, highnumsplits, -1);
}

void mesh::split(int n)
{
    errorifloaded();
    rawmeshptr->split(n);
}

void mesh::move(int physreg, expression u)
{
    errorifnotloaded();
    rawmeshptr->getphysicalregions()->errorundefined({physreg});
    rawmeshptr->gethadaptedpointer()->getphysicalregions()->errorundefined({physreg});
    
    if (rawmeshptr->getmeshnumber() == 0)
        rawmeshptr->move(physreg, u);
    rawmeshptr->gethadaptedpointer()->move(physreg, u);
}

void mesh::move(expression u)
{
    errorifnotloaded();
    if (rawmeshptr->getmeshnumber() == 0)
        rawmeshptr->move(-1, u);
    rawmeshptr->gethadaptedpointer()->move(-1, u);
}
        
void mesh::shift(int physreg, double x, double y, double z)
{
    errorifnotloaded();
    rawmeshptr->getphysicalregions()->errorundefined({physreg});
    rawmeshptr->gethadaptedpointer()->getphysicalregions()->errorundefined({physreg});
    
    rawmeshptr->shift(physreg, x, y, z);
    rawmeshptr->gethadaptedpointer()->shift(physreg, x, y, z);
}

void mesh::shift(double x, double y, double z)
{
    errorifnotloaded();
    rawmeshptr->shift(-1, x, y, z);
    rawmeshptr->gethadaptedpointer()->shift(-1, x, y, z);
}

void mesh::rotate(int physreg, double ax, double ay, double az)
{
    errorifnotloaded();
    rawmeshptr->getphysicalregions()->errorundefined({physreg});
    rawmeshptr->gethadaptedpointer()->getphysicalregions()->errorundefined({physreg});
    
    rawmeshptr->rotate(physreg, ax, ay, az);
    rawmeshptr->gethadaptedpointer()->rotate(physreg, ax, ay, az);
}

void mesh::rotate(double ax, double ay, double az)
{
    errorifnotloaded();
    rawmeshptr->rotate(-1, ax, ay, az);
    rawmeshptr->gethadaptedpointer()->rotate(-1, ax, ay, az);
}

void mesh::scale(int physreg, double x, double y, double z)
{
    errorifnotloaded();
    rawmeshptr->getphysicalregions()->errorundefined({physreg});
    rawmeshptr->gethadaptedpointer()->getphysicalregions()->errorundefined({physreg});
    
    rawmeshptr->scale(physreg, x, y, z);
    rawmeshptr->gethadaptedpointer()->scale(physreg, x, y, z);
}

void mesh::scale(double x, double y, double z)
{
    errorifnotloaded();
    rawmeshptr->scale(-1, x, y, z);
    rawmeshptr->gethadaptedpointer()->scale(-1, x, y, z);
}

int mesh::getmeshdimension(void)
{
    errorifnotloaded();
    return rawmeshptr->getmeshdimension();
}

std::vector<int> mesh::getphysicalregionnumbers(int dim)
{
    errorifnotloaded();
    return rawmeshptr->getphysicalregionnumbers(dim);
}

void mesh::selectskin(int newphysreg, int physregtoskin)
{
    if (physregtoskin < 0)
    {
        std::cout << "Error in 'mesh' object: expected a positive physical region number" << std::endl;
        abort();
    }

    errorifloaded();
    rawmeshptr->selectskin(newphysreg, physregtoskin);
}

void mesh::selectskin(int newphysreg)
{
    errorifloaded();
    rawmeshptr->selectskin(newphysreg, -1);
}

void mesh::selectbox(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit)
{
    if (physregtobox < 0)
    {
        std::cout << "Error in 'mesh' object: expected a positive physical region number" << std::endl;
        abort();
    }
    
    errorifloaded();
    rawmeshptr->selectbox(newphysreg, physregtobox, selecteddim, boxlimit);
}

void mesh::selectbox(int newphysreg, int selecteddim, std::vector<double> boxlimit)
{
    errorifloaded();
    rawmeshptr->selectbox(newphysreg, -1, selecteddim, boxlimit);
}

void mesh::selectsphere(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius)
{
    if (physregtosphere < 0)
    {
        std::cout << "Error in 'mesh' object: expected a positive physical region number" << std::endl;
        abort();
    }
    
    errorifloaded();
    rawmeshptr->selectsphere(newphysreg, physregtosphere, selecteddim, centercoords, radius);
}

void mesh::selectsphere(int newphysreg, int selecteddim, std::vector<double> centercoords, double radius)
{
    errorifloaded();
    rawmeshptr->selectsphere(newphysreg, -1, selecteddim, centercoords, radius);
}

void mesh::selectlayer(int newphysreg, int physregtoselectfrom, int physregtostartgrowth, int numlayers)
{
    if (physregtoselectfrom < 0)
    {
        std::cout << "Error in 'mesh' object: expected a positive physical region number" << std::endl;
        abort();
    }
    
    errorifloaded();
    rawmeshptr->selectlayer(newphysreg, physregtoselectfrom, physregtostartgrowth, numlayers);
}

void mesh::selectlayer(int newphysreg, int physregtostartgrowth, int numlayers)
{
    errorifloaded();
    rawmeshptr->selectlayer(newphysreg, -1, physregtostartgrowth, numlayers);
}

void mesh::selectexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude)
{
    if (physregtoexcludefrom < 0)
    {
        std::cout << "Error in 'mesh' object: expected a positive physical region number" << std::endl;
        abort();
    }
    
    errorifloaded();
    rawmeshptr->selectexclusion(newphysreg, physregtoexcludefrom, physregstoexclude);
}

void mesh::selectexclusion(int newphysreg, std::vector<int> physregstoexclude)
{
    errorifloaded();
    rawmeshptr->selectexclusion(newphysreg, -1, physregstoexclude);
}

void mesh::selectanynode(int newphysreg, int physregtoselectfrom)
{
    if (physregtoselectfrom < 0)
    {
        std::cout << "Error in 'mesh' object: expected a positive physical region number" << std::endl;
        abort();
    }
    
    errorifloaded();
    rawmeshptr->selectanynode(newphysreg, physregtoselectfrom);
}

void mesh::selectanynode(int newphysreg)
{
    errorifloaded();
    rawmeshptr->selectanynode(newphysreg, -1);
}

void mesh::use(void)
{
    errorifnotloaded();
    universe::mymesh = rawmeshptr->gethadaptedpointer();
}

