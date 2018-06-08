#include "mesh.h"


void mesh::readfromfile(std::string name)
{
	filename = name;

	if (name.length() >= 5 && name.compare(name.size()-4,4,".msh") == 0)
		gmshinterface::readfromfile(name, mynodes, myelements, myphysicalregions);
    else
    {
		std::cout << "Error: file '" << name << "' has either no extension or it is not supported." << std::endl << "Currently supported: GMSH .msh" << std::endl;
		abort();
	}
}

void mesh::writetofile(std::string name)
{
	if (name.length() >= 5 && name.compare(name.size()-4,4,".msh") == 0)
		gmshinterface::writetofile(name, mynodes, myelements, myphysicalregions, mydisjointregions);
    else
    {
		std::cout << "Error: file '" << name << "' has either no extension or it is not supported." << std::endl << "Currently supported: GMSH .msh" << std::endl;
		abort();
	}
}

void mesh::sortbybarycenters(void)
{
	for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        if (myelements.count(elementtypenumber) == 0)
            continue;
        // Sort the elements of the current type by their barycenters:
        std::vector<int> renumberingvector = myelements.sortbybarycenters(elementtypenumber);
        
        // Renumber the elements in the physical regions:
        for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
        {
            physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
            currentphysicalregion->renumberelements(elementtypenumber, renumberingvector);
        }
    }
}

void mesh::removeduplicates(void)
{
	for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        if (myelements.count(elementtypenumber) == 0)
            continue;
        // Remove the element duplicates:
        std::vector<int> renumberingvector = myelements.removeduplicates(elementtypenumber);
        
        // Renumber the elements in the physical regions:
        for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
        {
            physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
            currentphysicalregion->renumberelements(elementtypenumber, renumberingvector);
        }
    }
    // Remove the duplicated elements in every physical region:
    for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
        currentphysicalregion->removeduplicatedelements();
    }
}

void mesh::printcount(void)
{
	for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
	{
		if (myelements.count(elementtypenumber) == 0)
			continue;
			
		// Sparselizard does not support pyramids for the moment:
		if (elementtypenumber == 7)
		{
			std::cout << "Error in 'mesh' object: sparselizard does not support pyramid elements for the moment" << std::endl;
			abort();
		}
        
		element elementobject(elementtypenumber);
		std::string elementname = elementobject.gettypenameconjugation(myelements.count(elementtypenumber));
        // For nodes there is no curvature order:
        if (elementtypenumber == 0)
            std::cout << "Extracted " << myelements.count(elementtypenumber) << " " << "nodes" << std::endl;
        else
            std::cout << "Extracted " << myelements.count(elementtypenumber) << " " << elementname << " with curvature order " << myelements.getcurvatureorder() << std::endl;
	}
}

mesh::mesh(void) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions)
{
    universe::mymesh = this;
}

mesh::mesh(std::string filename, int verbosity) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions)
{
    universe::mymesh = this;
    load(filename, verbosity);
}

nodes* mesh::getnodes(void) {return &mynodes;}
elements* mesh::getelements(void) {return &myelements;}
physicalregions* mesh::getphysicalregions(void) {return &myphysicalregions;}
disjointregions* mesh::getdisjointregions(void) {return &mydisjointregions;}

void mesh::load(std::string name, int verbosity)
{
	///// Reset all memory of the mesh object:
	mynodes = nodes();
	mydisjointregions = disjointregions();
	myphysicalregions = physicalregions(mydisjointregions);
 	myelements = elements(mynodes, myphysicalregions, mydisjointregions);
	///// Memory is reset


    if (verbosity > 0)
        std::cout << "Loading mesh from file '" << name << "'" << std::endl;
    
	wallclock loadtime;
    
    readfromfile(name);
	myelements.explode();
	sortbybarycenters();
	removeduplicates();
    
    myelements.definedisjointregions();
    // The reordering is stable and the elements are thus still ordered 
    // by barycenter coordinates in every disjoint region!
    myelements.reorderbydisjointregions();
    myelements.definedisjointregionsranges();
    
    // Define the physical regions based on the disjoint regions they contain:
    for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
        currentphysicalregion->definewithdisjointregions();
    }
    
    myelements.tostandardorientation();
	myelements.orient();
    
    if (verbosity > 0)
    {
        printcount();
        loadtime.print("Time to read the mesh: ");
    }
}

void mesh::write(std::string name, int verbosity)
{
    if (verbosity > 0)
        std::cout << "Writing mesh to file '" << name << "'" << std::endl;
    
	wallclock writetime;
	
    writetofile(name);

    if (verbosity > 0)
    {
        std::cout << "Time to write the mesh: ";
        writetime.print();
    }
}

void mesh::shift(double x, double y, double z)
{
    mynodes.shift(x, y, z);
}

void mesh::rotate(double ax, double ay, double az)
{
    mynodes.rotate(ax, ay, az);
}

int mesh::getmeshdimension(void)
{
	int maxelementdimension = -1;

    // Iterate on all element types in 'myelements':
	for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
	{
		// If there is no element of that type we skip it:
		if (myelements.count(elementtypenumber) == 0)
			continue;
		
		element elementobject(elementtypenumber);
		if (elementobject.getelementdimension() > maxelementdimension)
			maxelementdimension = elementobject.getelementdimension();
	}
	return maxelementdimension;
}


void mesh::writewithdisjointregions(std::string name)
{
    // Create a new 'physicalregions' object where the ith 
    // physical region contains only disjoint region i. 
    physicalregions tempphysicalregions(mydisjointregions);
    for (int disjreg = 0; disjreg < mydisjointregions.count(); disjreg++)
    {
        physicalregion* currentphysicalregion = tempphysicalregions.get(disjreg+1); // +1 because GMSH does not accept 0 as physical region number
        currentphysicalregion->setdisjointregions({disjreg});
    }
    // Save to GMSH format:
    gmshinterface::writetofile(name, mynodes, myelements, tempphysicalregions, mydisjointregions);
}

void mesh::printphysicalregions(void)
{
    std::cout << std::endl << "Element dimension and disjoint region list in every physical region:" << std::endl;
    for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
        std::cout << myphysicalregions.getnumber(physregindex) << " (" << currentphysicalregion->getelementdimension() << "D) : ";
        
        std::vector<int> disjointregionlist = currentphysicalregion->getdisjointregions();
        // Print the disjoint region list:
        for (int i = 0; i < disjointregionlist.size(); i++)
            std::cout << disjointregionlist[i] << " ";
        
        std::cout << std::endl;
    }
}

void mesh::printdisjointregions(void)
{
    std::cout << std::endl << "Element type, #elements and physical region list in every disjoint region:" << std::endl;
    for (int disjreg = 0; disjreg < mydisjointregions.count(); disjreg++)
    {
        std::cout << disjreg << " (" << mydisjointregions.getelementtypenumber(disjreg) << ", " << mydisjointregions.countelements(disjreg) << ") : ";
        
        for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
        {
            if (mydisjointregions.isinphysicalregion(disjreg, physregindex))
                std::cout << myphysicalregions.getnumber(physregindex) << " ";
        }
        std::cout << std::endl;
    }
}






