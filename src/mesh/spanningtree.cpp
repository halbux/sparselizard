#include "spanningtree.h"
#include "universe.h"
#include "slexceptions.h"


spanningtree::spanningtree(std::vector<int> physregs)
{
    if (universe::myrawmesh == NULL)
    {
        throw slexception( "Error in 'spanningtree' object: cannot define a spanning tree before the mesh is loaded" );
    }
    
    universe::getrawmesh()->getphysicalregions()->errorundefined(physregs);
    
    myrawspantree = std::shared_ptr<rawspanningtree>(new rawspanningtree(physregs));
}

int spanningtree::countedgesintree(void)
{
    return myrawspantree->countedgesintree();
}

std::vector<int> spanningtree::getedgesintree(void)
{
    return myrawspantree->getedgesintree();
}

std::shared_ptr<rawspanningtree> spanningtree::getpointer(void)
{
    return myrawspantree;
}

void spanningtree::write(std::string filename)
{
    myrawspantree->write(filename);
}

