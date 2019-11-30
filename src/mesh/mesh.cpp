#include "mesh.h"


void mesh::readfromfile(std::string name)
{
    if (name.length() >= 5 && name.compare(name.size()-4,4,".msh") == 0)
    {
        gmshinterface::readfromfile(name, mynodes, myelements, myphysicalregions);
        return;    
    }
    if (name.length() >= 5 && name.compare(name.size()-4,4,".nas") == 0)
    {
        nastraninterface::readfromfile(name, mynodes, myelements, myphysicalregions);
        return;    
    }

    std::cout << "Error: file '" << name << "' cannot be read by the legacy mesh reader." << std::endl;
    std::cout << "Use the petsc mesh reader instead or use the GMSH .msh or Nastran .nas format." << std::endl;
    abort();
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

mesh::mesh(void) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions), myregiondefiner(mynodes, myelements, myphysicalregions)
{
    universe::mymesh = this;
}

mesh::mesh(std::string filename, int verbosity, bool legacyreader) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions), myregiondefiner(mynodes, myelements, myphysicalregions)
{
    universe::mymesh = this;
    load(filename, verbosity, legacyreader);
}

mesh::mesh(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions), myregiondefiner(mynodes, myelements, myphysicalregions)
{
    universe::mymesh = this;
    load(mergeduplicates, meshfiles, verbosity);
}

mesh::mesh(std::vector<shape> inputshapes, int verbosity) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions), myregiondefiner(mynodes, myelements, myphysicalregions)
{
    universe::mymesh = this;
    load(inputshapes, verbosity);
}

nodes* mesh::getnodes(void) {return &mynodes;}
elements* mesh::getelements(void) {return &myelements;}
physicalregions* mesh::getphysicalregions(void) {return &myphysicalregions;}
disjointregions* mesh::getdisjointregions(void) {return &mydisjointregions;}

void mesh::load(std::string name, int verbosity, bool legacyreader)
{
    ///// Reset all memory of the mesh object:
    mynodes = nodes();
    mydisjointregions = disjointregions();
    myphysicalregions = physicalregions(mydisjointregions);
    myelements = elements(mynodes, myphysicalregions, mydisjointregions);
    ///// Memory is reset


    filename = name;

    if (verbosity > 0)
        std::cout << "Loading mesh from file '" << name << "'" << std::endl;
    
    wallclock loadtime;
    
    if (legacyreader)
        readfromfile(name);
    else
    {
        petscmesh pmesh(name);
        pmesh.extract(mynodes, myelements, myphysicalregions);
    }
    myelements.explode();
    sortbybarycenters();
    removeduplicates();
    myregiondefiner.defineregions();
    
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
        printcount();
    if (verbosity > 1)
        printelementsinphysicalregions();
    if (verbosity > 0)
        loadtime.print("Time to load the mesh: ");

    // Make sure axisymmetry is valid for this mesh:    
    if (universe::isaxisymmetric && getmeshdimension() != 2)
    {
        std::cout << "Error in 'mesh' object: axisymmetry is only allowed for 2D problems" << std::endl;
        abort();
    }
}

void mesh::load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
{
    int numfiles = meshfiles.size();
    if (numfiles == 0)
    {
        std::cout << "Error in 'mesh' object: expected at least one mesh file to load" << std::endl;
        abort();
    }
    
    std::vector< std::vector<std::vector<shape>> > allshapes(numfiles, std::vector<std::vector<shape>>(0));
    
    int meshdim = -1;

    int maxphysreg = 0;
    for (int i = 0; i < numfiles; i++)
    {
        allshapes[i] = mathop::loadshape(meshfiles[i]);
        
        for (int d = 0; d < 4; d++)
        {
            int curcount = allshapes[i][d].size();
        
            if (i == 0 && curcount > 0)
                meshdim = d;
                
            if (i > 0 && allshapes[i][meshdim].size() == 0 || i > 0 && d > meshdim && allshapes[i][d].size() > 0)
            {
                std::cout << "Error in 'mesh' object: expected a " << meshdim << "D mesh in '" << meshfiles[i] << "'" << std::endl;
                abort();
            }
        
            for (int s = 0; s < curcount; s++)
            {
                int curphysreg = allshapes[i][d][s].getphysicalregion();
                if (maxphysreg < curphysreg)
                    maxphysreg = curphysreg;
            }
        }
    }
    
    // Amount to shift each mesh file to make sure they do not overlap:
    std::vector<double> shiftvec(numfiles,0.0);
    
    std::vector<shape> flat = {};
    std::vector<double> prevcoordbounds;
    for (int i = 0; i < numfiles; i++)
    {
        shape curunion("union", maxphysreg+1+i, allshapes[i][meshdim]);
        
        if (mergeduplicates == false)
        {
            std::vector<double> coords = curunion.getcoords();
            std::vector<double> curcoordbounds = myalgorithm::getcoordbounds(coords);
            if (i > 0)
                shiftvec[i] = shiftvec[i-1] + prevcoordbounds[1]-curcoordbounds[0] + 0.1*(prevcoordbounds[1]-prevcoordbounds[0]);
            prevcoordbounds = curcoordbounds;
        }
        
        // Do not shift recursively.
        if (mergeduplicates == false)
            curunion.getpointer()->shift(shiftvec[i],0,0, false);
        flat.push_back(curunion);
    
        for (int d = 0; d < 4; d++)
        {
            for (int s = 0; s < allshapes[i][d].size(); s++)
            {
                if (mergeduplicates == false)
                    allshapes[i][d][s].getpointer()->shift(shiftvec[i],0,0, false);
                flat.push_back(allshapes[i][d][s]);
            }
        }
    }
    
    this->load(flat, verbosity);
    
    // Shift back to the original position:
    if (mergeduplicates == false)
    {
        for (int i = 0; i < numfiles; i++)
            this->shift(maxphysreg+1+i, -shiftvec[i],0,0);
    }
}

void mesh::load(std::vector<shape> inputshapes, int verbosity)
{
    ///// Reset all memory of the mesh object:
    mynodes = nodes();
    mydisjointregions = disjointregions();
    myphysicalregions = physicalregions(mydisjointregions);
    myelements = elements(mynodes, myphysicalregions, mydisjointregions);
    filename = "";
     ///// Memory is reset


    if (verbosity > 0)
        std::cout << "Loading mesh from " << inputshapes.size() << " shapes" << std::endl;
    
    wallclock loadtime;
    

    ///// Transfer the mesh from every shape to the corresponding objects:
    // Curvature order for shapes is 1 for now:
    int curvatureorder = 1;

    if (inputshapes.size() == 0)
    {
        std::cout << "Error in 'mesh' object: provided an empty vector of shapes" << std::endl;
        abort();
    }

    // Make sure all shapes have a valid physical region number:
    for (int i = 0; i < inputshapes.size(); i++)
    {
        if (inputshapes[i].getphysicalregion() <= 0)
        {
            std::cout << "Error in 'mesh' object: physical region was undefined (negative or zero) for shape '" << inputshapes[i].getname() << "'" << std::endl;
            abort();
        }
    }

    // Get the number of nodes for preallocation:
    int numberofnodes = 0;
    for (int i = 0; i < inputshapes.size(); i++)
        numberofnodes += ( inputshapes[i].getpointer()->getcoords() )->size()/3;
    mynodes.setnumber(numberofnodes);
    std::vector<double>* nodecoordinates = mynodes.getcoordinates();


    // The node numbers must be shifted from a shape to the other to avoid same numbers:
    int offset = 0;
    // Loop on every input shape:
    for (int i = 0; i < inputshapes.size(); i++)
    {
        // Append the nodes:
        std::vector<double>* nodecoords = inputshapes[i].getpointer()->getcoords();
        for (int j = 0; j < nodecoords->size(); j++)
            nodecoordinates->at(3*offset+j) = nodecoords->at(j);

        // Append the elements:
        int physreg = inputshapes[i].getpointer()->getphysicalregion();
        physicalregion* currentphysicalregion = myphysicalregions.get(physreg);
        std::vector<std::vector<int>>* elems = inputshapes[i].getpointer()->getelems();
        // Loop on all element types:
        for (int typenum = 0; typenum < elems->size(); typenum++)
        {
            element currentelem(typenum, curvatureorder);
            // Number of nodes in every element of current type number:
            int numnodesinelem = currentelem.countcurvednodes();

            std::vector<int> nodesincurrentelement(numnodesinelem);

            for (int e = 0; e < (elems->at(typenum)).size()/numnodesinelem; e++)
            {
                for (int m = 0; m < numnodesinelem; m++)
                    nodesincurrentelement[m] = (elems->at(typenum))[e*numnodesinelem+m] + offset;

                int elementindexincurrenttype = myelements.add(typenum, curvatureorder, nodesincurrentelement);
                currentphysicalregion->addelement(typenum, elementindexincurrenttype);
            }
        }

        offset += nodecoords->size()/3;
    }
    ///// Mesh is transferred


    myelements.explode();
    sortbybarycenters();
    removeduplicates();
    myregiondefiner.defineregions();
    
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
        printcount();
    if (verbosity > 1)
        printelementsinphysicalregions();
    if (verbosity > 0)
        loadtime.print("Time to load the mesh: ");

    // Make sure axisymmetry is valid for this mesh:    
    if (universe::isaxisymmetric && getmeshdimension() != 2)
    {
        std::cout << "Error in 'mesh' object: axisymmetry is only allowed for 2D problems" << std::endl;
        abort();
    }
}


void mesh::write(std::string name, int verbosity)
{
    if (verbosity > 0)
        std::cout << "Writing mesh to file '" << name << "'" << std::endl;
    
    wallclock writetime;
    
    writetofile(name);

    if (verbosity > 0)
        writetime.print("Time to write the mesh: ");
}

void mesh::shift(int physreg, double x, double y, double z)
{
    mynodes.shift(physreg, x, y, z);
    myelements.shift(physreg, x, y, z);
}

void mesh::shift(double x, double y, double z)
{
    mynodes.shift(-1, x, y, z);
    myelements.shift(-1, x, y, z);
}

void mesh::rotate(int physreg, double ax, double ay, double az)
{
    mynodes.rotate(physreg, ax, ay, az);
    myelements.rotate(physreg, ax, ay, az);
}

void mesh::rotate(double ax, double ay, double az)
{
    mynodes.rotate(-1, ax, ay, az);
    myelements.rotate(-1, ax, ay, az);
}

void mesh::scale(int physreg, double x, double y, double z)
{
    mynodes.scale(physreg, x, y, z);
    myelements.scale(physreg, x, y, z);
}

void mesh::scale(double x, double y, double z)
{
    mynodes.scale(-1, x, y, z);
    myelements.scale(-1, x, y, z);
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

void mesh::regionskin(int newphysreg, int physregtoskin)
{
    myregiondefiner.regionskin(newphysreg, physregtoskin);
}

void mesh::boxselection(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit)
{
    myregiondefiner.boxselection(newphysreg, selecteddim, boxlimit, physregtobox);
}

void mesh::sphereselection(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius)
{
    myregiondefiner.sphereselection(newphysreg, selecteddim, centercoords, radius, physregtosphere);
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

void mesh::printelementsinphysicalregions(bool isdebug)
{
    int numphysregs = myphysicalregions.count();
    
    std::cout << "Extracted " << numphysregs << " physical region"+myalgorithm::getplurals(numphysregs)+":" << std::endl;
    for (int physregindex = 0; physregindex < numphysregs; physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
        
        if (isdebug == false)
        {
            int numelems = currentphysicalregion->countelements();
            std::cout << myphysicalregions.getnumber(physregindex) << " (" << numelems << " " << currentphysicalregion->getelementdimension() << "D element"+myalgorithm::getplurals(numelems)+")";
        }
        else
        {
            std::cout << myphysicalregions.getnumber(physregindex) << " (" << currentphysicalregion->getelementdimension() << "D)";
            
            std::vector<std::vector<int>>* elemlist = currentphysicalregion->getelementlist();

            for (int i = 0; i < elemlist->size(); i++)
            {
                if (elemlist->at(i).size() > 0)
                {
                    std::cout << " (" << i << "): ";
                    for (int j = 0; j < elemlist->at(i).size(); j++)
                        std::cout << elemlist->at(i)[j] << " ";
                }
            }
        }

        std::cout << std::endl;
    }
}

