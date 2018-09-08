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
	isedgeintree = std::vector<bool>(myelements->count(1),false);


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

				isedgeintree[edgesinsubtree[currentsubtree][i]] = true;
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

		// Get the disjoint edge region number of the current candidate:
        int currentdisjreg = myelements->getdisjointregion(1, currentedge);

		if (isedgeintree[currentedge])
			continue;

		// Get the other node in the current edge:
		int nextnode = myelements->getsubelement(0, 1, currentedge, 0);
		if (nodenumber == nextnode)
			nextnode = myelements->getsubelement(0, 1, currentedge, 1);

		// Add the edge to the tree if it does not create a loop:
		if (isnodeintree[nextnode] == false)
		{
			isnodeintree[nextnode] = true;
			isedgeintree[currentedge] = true;
			numberofedgesintree++;
			
			growtree(nextnode);
		}
	}
}

spanningtree::spanningtree(std::vector<int> physregs)
{
	myelements = universe::mymesh->getelements();
	mydisjointregions = universe::mymesh->getdisjointregions();


	// Get a vector with all disjoint edge regions in the physical regions provided:
	isprioritydisjointregion = std::vector<bool>(mydisjointregions->count(), false);

	for (int i = 0; i < physregs.size(); i++)
	{
		// Get all disjoint edge regions in the current physical region:
		std::vector<int> edgedisjregs = ((universe::mymesh->getphysicalregions())->get(physregs[i]))->getdisjointregions(1);

		for (int j = 0; j < edgedisjregs.size(); j++)
			isprioritydisjointregion[edgedisjregs[j]] = true;
	}

	// Grow all subtrees (they share no node and no edge):
	growsubtrees();
	// Connect the subtrees together:
	connectsubtrees();
}

bool spanningtree::isintree(int index, int disjreg)
{
	bool output = false;
	if (mydisjointregions->getelementtypenumber(disjreg) == 1)
		output = isedgeintree[mydisjointregions->getrangebegin(disjreg)+index];

	return output;
}

int spanningtree::countedgesintree(int disjreg)
{
	int output = 0;
	if (mydisjointregions->getelementtypenumber(disjreg) == 1)
	{
		int derstartindex = mydisjointregions->getrangebegin(disjreg);
		int derendindex = mydisjointregions->getrangeend(disjreg);

		for (int i = derstartindex; i <= derendindex; i++)
		{
			if (isedgeintree[i])
				output++;
		}
	}

	return output;
}

int spanningtree::countedgesintree(void)
{
	return numberofedgesintree;
}

std::vector<int> spanningtree::getedgesintree(void)
{
	std::vector<int> output(numberofedgesintree);

	int index = 0;
	for (int i = 0; i < isedgeintree.size(); i++)
	{
		if (isedgeintree[i])
		{
			output[index] = i;
			index++;
		}
	}

	return output;
}

void spanningtree::write(std::string filename)
{
    // Make sure the filename includes the extension:
    if (filename.size() < 5 || filename.substr(filename.size()-4,4) != ".pos")
    {
        std::cout << "Error in 'spanningtree' object: cannot write to file '" << filename << "' (unknown or missing file extension)" << std::endl;
        abort();
    }

	nodes* mynodes = universe::mymesh->getnodes();
	std::vector<double>* nodecoords = mynodes->getcoordinates();
	
	densematrix xcoords(numberofedgesintree,2), ycoords(numberofedgesintree,2), zcoords(numberofedgesintree,2);

	double* xvals = xcoords.getvalues();
	double* yvals = ycoords.getvalues();	
	double* zvals = zcoords.getvalues();		

	int index = 0;
	for (int i = 0; i < isedgeintree.size(); i++)
	{
		if (isedgeintree[i])
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

    // Write the header:
    gmshinterface::openview(filename, filename, 0, true);
    // Write the data:
	gmshinterface::appendtoview(filename, 1, xcoords, ycoords, zcoords, densematrix(numberofedgesintree,2,1.0));
    // Close view:
    gmshinterface::closeview(filename);

}




