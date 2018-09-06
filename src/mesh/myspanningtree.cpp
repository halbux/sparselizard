#include "myspanningtree.h"


void myspanningtree::growsubtrees(void)
{
	// Initialise:
	insubtree = std::vector<int>(myelements->count(1), -1);
	isnodeintree = std::vector<bool>(myelements->count(0),false);

	numberofsubtrees = 0;
	for (int i = 0; i < ordereddisjointedgeregions.size(); i++)
	{
		int currentder = ordereddisjointedgeregions[i];
		int edgerangebegin = mydisjointregions->getrangebegin(currentder);
		for (int j = 0; j < mydisjointregions->countelements(currentder); j++)
		{
			int startnode = myelements->getsubelement(0, 1, edgerangebegin+j, 0);
			int endnode = myelements->getsubelement(0, 1, edgerangebegin+j, 1);

			// If at least the current edge can be added to the subtree:
			if (isnodeintree[startnode] == false && isnodeintree[endnode] == false)
			{
				growsubtree(currentder, startnode, numberofsubtrees);
				numberofsubtrees++;
			}
		}
	}
}

void myspanningtree::growsubtree(int edgedisjreg, int nodenumber, int subtreenumber)
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

		// Skip edge if not in the requested disjoint region:
		if (currentder != edgedisjreg)
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
			growsubtree(edgedisjreg, nextnode, subtreenumber);
		}
	}
}

void myspanningtree::connectsubtrees(void)
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


void myspanningtree::growtree(int nodenumber)
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

myspanningtree::myspanningtree(void)
{
	myelements = universe::mymesh->getelements();
	mydisjointregions = universe::mymesh->getdisjointregions();


	// Get all disjoint edge regions in the mesh:
	std::vector<int> alledgedisjregs = mydisjointregions->get(1);
	int numedgedisjregs = alledgedisjregs.size();


	// Create a vector with the dimension of the geometrical shape in which every disjoint edge region is (1D, 2D, 3D).
	// For that take the intersection of all physical regions defining a disjoint edge region.
	std::vector<int> derdims(numedgedisjregs, 0);

	for (int e = 0; e < numedgedisjregs; e++)
	{
		// Get the physical regions defining the current disjoint edge region:
		std::vector<int> physregsdefiningder = mydisjointregions->getphysicalregions(alledgedisjregs[e]);

		// Take the intersection of all disjoint regions in the physregs defining the current disjoint edge region:
		std::vector<int> disjregs = {};
		for (int i = 0; i < physregsdefiningder.size(); i++)
		{
		    // Get all disjoint regions with -1:
		    std::vector<int> disjregsinthisphysreg = universe::mymesh->getphysicalregions()->get(physregsdefiningder[i])->getdisjointregions(-1);
		    
		    if (i > 0)
		        disjregs = myalgorithm::intersect(disjregs, disjregsinthisphysreg);
		    else
		        disjregs = disjregsinthisphysreg;
		}

		// The dim we are looking for is the max dim of any disjoint region in the intersected list:
		for (int i = 0; i < disjregs.size(); i++)
		{
			if (derdims[e] < mydisjointregions->getelementdimension(disjregs[i]))
				derdims[e] = mydisjointregions->getelementdimension(disjregs[i]);
		}
    }
	

	// Sort the disjoint edge regions with increasing dimension (can only be 1D, 2D or 3D):
	ordereddisjointedgeregions = std::vector<int>(numedgedisjregs);

	int index = 0;
	for (int d = 1; d <= 3; d++)
	{
		for (int i = 0; i < numedgedisjregs; i++)
		{
			if (derdims[i] == d)
			{
				ordereddisjointedgeregions[index] = alledgedisjregs[i];
				index++;
			}
		}
	}


	// Grow all subtrees (they share no node and no edge):
	growsubtrees();
	// Connect the subtrees together:
	connectsubtrees();
}

bool myspanningtree::isintree(int index, int disjreg)
{
	bool output = false;
	if (mydisjointregions->getelementtypenumber(disjreg) == 1)
		output = isedgeintree[mydisjointregions->getrangebegin(disjreg)+index];

	return output;
}

int myspanningtree::countedgesintree(int disjreg)
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

int myspanningtree::countedgesintree(void)
{
	return numberofedgesintree;
}

std::vector<int> myspanningtree::getedgesintree(void)
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

void myspanningtree::write(std::string filename)
{
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




