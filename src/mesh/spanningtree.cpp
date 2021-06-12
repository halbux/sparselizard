#include "spanningtree.h"


void spanningtree::growsubtrees(void)
{
    // Initialise:
    insubtree = std::vector<int>(myelements->count(1), -1);
    isnodeintree = std::vector<bool>(myelements->count(0),false);

    numberofsubtrees = 0;
    for (int i = 0; i < isprioritydisjointregion.size(); i++)
    {
        if (isprioritydisjointregion[i] == false)
            continue;

        int edgerangebegin = mydisjointregions->getrangebegin(i);
        for (int j = 0; j < mydisjointregions->countelements(i); j++)
        {
            int startnode = myelements->getsubelement(0, 1, edgerangebegin+j, 0);
            int endnode = myelements->getsubelement(0, 1, edgerangebegin+j, 1);

            // If at least the current edge can be added to the subtree:
            if (isnodeintree[startnode] == false && isnodeintree[endnode] == false)
            {
                growsubtree(startnode, numberofsubtrees);
                numberofsubtrees++;
            }
        }
    }
}

void spanningtree::growsubtree(int nodenumber, int subtreenumber)
{
    isnodeintree[nodenumber] = true;

    // Get the list of edges touching the node 'nodenumber':
    std::vector<int> edgecandidates = myelements->getedgesonnode(nodenumber);    

    // Loop on all edge candidates:
    for (int i = 0; i < edgecandidates.size(); i++)
    {
        int currentedge = edgecandidates[i];

        // Get the disjoint edge region number of the current candidate:
        int currentder = myelements->getdisjointregion(1, currentedge);

        // Skip edge if not in the priority disjoint edge regions:
        if (isprioritydisjointregion[currentder] == false)
            continue;

        // Get the other node in the current edge:
        int nextnode = myelements->getsubelement(0, 1, currentedge, 0);
        if (nodenumber == nextnode)
            nextnode = myelements->getsubelement(0, 1, currentedge, 1);

        // Add the edge to the tree if it does not create a loop:
        if (isnodeintree[nextnode] == false)
        {
            isnodeintree[nextnode] = true;
            insubtree[currentedge] = subtreenumber;
            growsubtree(nextnode, subtreenumber);
        }
    }
}

void spanningtree::connectsubtrees(void)
{
    isnodeintree = std::vector<bool>(myelements->count(0),false);

    isedgeintree = std::shared_ptr<bool>(new bool[myelements->count(1)]);
    isedgeintreeptr = isedgeintree.get();
    // Initialise to all false:
    for (int i = 0; i < myelements->count(1); i++)
        isedgeintreeptr[i] = false;


    ///// Group the edges of each subtree in a vector:


    // Count the number of edges in each subtree for preallocation:
    std::vector<int> numedgesinsubtree(numberofsubtrees,0);
    for (int i = 0; i < insubtree.size(); i++)
    {
        if (insubtree[i] != -1)
            numedgesinsubtree[insubtree[i]]++;
    }
    // Preallocate:
    edgesinsubtree = std::vector<std::vector<int>>(numberofsubtrees);
    for (int i = 0; i < numberofsubtrees; i++)
        edgesinsubtree[i].resize(numedgesinsubtree[i]);
    // Populate:
    std::vector<int> indexinsubtree(numberofsubtrees,0);
    for (int i = 0; i < insubtree.size(); i++)
    {
        if (insubtree[i] != -1)
        {
            edgesinsubtree[insubtree[i]][indexinsubtree[insubtree[i]]] = i;
            indexinsubtree[insubtree[i]]++;
        }
    }


    ///// Create the overall tree by connecting the subtrees:

    issubtreeintree = std::vector<bool>(numberofsubtrees,false);
    for (int i = 0; i < isnodeintree.size(); i++)
    {
        if (isnodeintree[i] == false)
            growtree(i);
    }

    
    ///// Clean all temporary containers:
    
    numberofsubtrees = 0;
    insubtree = {};
    edgesinsubtree = {};
    issubtreeintree = {};
    isnodeintree = {};

}


void spanningtree::growtree(int nodenumber)
{
    isnodeintree[nodenumber] = true;

    // Get the list of edges touching the node 'nodenumber':
    std::vector<int> edgecandidates = myelements->getedgesonnode(nodenumber);    

    // Add to the tree all not yet added subtrees on the edges candidates:
    for (int i = 0; i < edgecandidates.size(); i++)
    {
        int currentsubtree = insubtree[edgecandidates[i]];

        // Add the subtree to the tree if not already added:
        if (currentsubtree != -1 && issubtreeintree[currentsubtree] == false)
        {
            for (int i = 0; i < edgesinsubtree[currentsubtree].size(); i++)
            {
                int node1 = myelements->getsubelement(0, 1, edgesinsubtree[currentsubtree][i], 0);
                int node2 = myelements->getsubelement(0, 1, edgesinsubtree[currentsubtree][i], 1);

                isnodeintree[node1] = true;
                isnodeintree[node2] = true;

                isedgeintreeptr[edgesinsubtree[currentsubtree][i]] = true;
                numberofedgesintree++;
            }

            issubtreeintree[currentsubtree] = true;

            // Grow tree from every added edge:
            for (int i = 0; i < edgesinsubtree[currentsubtree].size(); i++)
            {
                int node1 = myelements->getsubelement(0, 1, edgesinsubtree[currentsubtree][i], 0);
                int node2 = myelements->getsubelement(0, 1, edgesinsubtree[currentsubtree][i], 1);

                growtree(node1);
                growtree(node2);
            }
        }
    }

    // Add to the tree all edges that do not create a loop:
    for (int i = 0; i < edgecandidates.size(); i++)
    {
        int currentedge = edgecandidates[i];

        if (isedgeintreeptr[currentedge])
            continue;

        // Get the other node in the current edge:
        int nextnode = myelements->getsubelement(0, 1, currentedge, 0);
        if (nodenumber == nextnode)
            nextnode = myelements->getsubelement(0, 1, currentedge, 1);

        // Add the edge to the tree if it does not create a loop:
        if (isnodeintree[nextnode] == false)
        {
            isnodeintree[nextnode] = true;
            isedgeintreeptr[currentedge] = true;
            numberofedgesintree++;
            
            growtree(nextnode);
        }
    }
}

void spanningtree::grow(void)
{
    // Get a vector with all disjoint edge regions in the physical regions provided:
    isprioritydisjointregion = std::vector<bool>(mydisjointregions->count(), false);

    for (int i = 0; i < startphysregs.size(); i++)
    {
        // Get all disjoint edge regions in the current physical region:
        std::vector<int> edgedisjregs = ((universe::mymesh->getphysicalregions())->get(startphysregs[i]))->getdisjointregions(1);

        for (int j = 0; j < edgedisjregs.size(); j++)
            isprioritydisjointregion[edgedisjregs[j]] = true;
    }

    // Grow all subtrees (they share no node and no edge):
    growsubtrees();
    // Connect the subtrees together:
    connectsubtrees();
}

void spanningtree::synchronize(void)
{
    if (issynchronizing || universe::mymesh->getmeshnumber() == mymeshnumber)
        return;
    issynchronizing = true;    


    myelements = universe::mymesh->getelements();
    mydisjointregions = universe::mymesh->getdisjointregions();
    
    // Reset structure:
    isedgeintree = NULL;
    isedgeintreeptr = NULL;
    numberofedgesintree = 0;
    isprioritydisjointregion = {};
    numberofsubtrees = 0;
    insubtree = {};
    edgesinsubtree = {};
    issubtreeintree = {};
    isnodeintree = {};
    
    grow();
    
    
    mymeshnumber = universe::mymesh->getmeshnumber();
    issynchronizing = false;
}

spanningtree::spanningtree(std::vector<int> physregs)
{
    universe::mymesh->getphysicalregions()->errorundefined(physregs);
    
    startphysregs = physregs;
    
    myelements = universe::mymesh->getelements();
    mydisjointregions = universe::mymesh->getdisjointregions();
    
    grow();
    
    mymeshnumber = universe::mymesh->getmeshnumber();
}

bool spanningtree::isintree(int index, int disjreg)
{
    synchronize();
    
    bool output = false;
    if (mydisjointregions->getelementtypenumber(disjreg) == 1)
        output = isedgeintreeptr[mydisjointregions->getrangebegin(disjreg)+index];

    return output;
}

int spanningtree::countedgesintree(int disjreg)
{
    synchronize();
    
    int output = 0;
    if (mydisjointregions->getelementtypenumber(disjreg) == 1)
    {
        int derstartindex = mydisjointregions->getrangebegin(disjreg);
        int derendindex = mydisjointregions->getrangeend(disjreg);

        for (int i = derstartindex; i <= derendindex; i++)
        {
            if (isedgeintreeptr[i])
                output++;
        }
    }

    return output;
}

int spanningtree::countedgesintree(void)
{
    synchronize();
    
    return numberofedgesintree;
}

std::vector<int> spanningtree::getedgesintree(void)
{
    synchronize();
    
    std::vector<int> output(numberofedgesintree);

    int index = 0;
    for (int i = 0; i < myelements->count(1); i++)
    {
        if (isedgeintreeptr[i])
        {
            output[index] = i;
            index++;
        }
    }

    return output;
}

void spanningtree::write(std::string filename)
{
    synchronize();
    
    nodes* mynodes = universe::mymesh->getnodes();
    std::vector<double>* nodecoords = mynodes->getcoordinates();
    
    densematrix xcoords(numberofedgesintree,2), ycoords(numberofedgesintree,2), zcoords(numberofedgesintree,2);

    double* xvals = xcoords.getvalues();
    double* yvals = ycoords.getvalues();    
    double* zvals = zcoords.getvalues();        

    int index = 0;
    for (int i = 0; i < myelements->count(1); i++)
    {
        if (isedgeintreeptr[i])
        {
            int node1 = myelements->getsubelement(0, 1, i, 0);
            int node2 = myelements->getsubelement(0, 1, i, 1);

            xvals[2*index+0] = nodecoords->at(3*node1+0);
            xvals[2*index+1] = nodecoords->at(3*node2+0);
            yvals[2*index+0] = nodecoords->at(3*node1+1);
            yvals[2*index+1] = nodecoords->at(3*node2+1);
            zvals[2*index+0] = nodecoords->at(3*node1+2);
            zvals[2*index+1] = nodecoords->at(3*node2+2);

            index++;
        }
    }
    
    // Write to file:
    iodata datatowrite(1, 1, true, {});
    datatowrite.addcoordinates(1, xcoords, ycoords, zcoords);
    datatowrite.adddata(1, {densematrix(numberofedgesintree, 2, 1.0)});
    
    iointerface::writetofile(filename, datatowrite);
}




