#include "rawmesh.h"
#include "geotools.h"


void rawmesh::splitmesh(void)
{
    if (mynumsplitrequested == 0)
        return;

    // Get all physical region numbers:
    std::vector<int> prn = myphysicalregions.getallnumbers();

    // Count the number of nodes per element:
    int co = myelements.getcurvatureorder();
    std::vector<int> ncn(8);
    for (int i = 0; i < 8; i++)
    {
        element myelem(i,co);
        ncn[i] = myelem.countcurvednodes();
    }
    
    // Loop on all physical regions:
    int numnodes = 0;
    std::vector< std::vector<std::vector<double>> > splitcoords(prn.size(), std::vector<std::vector<double>>(8, std::vector<double>(0)));
    for (int p = 0; p < prn.size(); p++)
    {
        physicalregion* curpr = myphysicalregions.getatindex(p);
        std::vector<std::vector<int>>* curelemlist = curpr->getelementlist();
    
        // Calculate the number of elements after split in the current region:
        std::vector<int> numsplit(8,0);
        for (int i = 0; i < 8; i++)
        {
            element myel(i,co);
            std::vector<int> fsc = myel.fullsplitcount(mynumsplitrequested);
        
            int ne = curelemlist->at(i).size();
            for (int j = 0; j < 8; j++)
                numsplit[j] += ne*fsc[j];
        }
        // Preallocate:
        for (int i = 0; i < 8; i++)
        {
            splitcoords[p][i] = std::vector<double>(3*ncn[i]*numsplit[i]);
            numnodes += ncn[i]*numsplit[i]; 
        }
                
        // Split each element type and group:
        std::vector<int> indexes(8,0);
        for (int i = 0; i < 8; i++)
        {
            int ne = curelemlist->at(i).size();
            if (ne == 0)
                continue;
            std::vector<double> coords(3*ncn[i]*ne);
            for (int e = 0; e < ne; e++)
            {
                int el = curelemlist->at(i)[e];
                std::vector<double> curcoords = myelements.getnodecoordinates(i,el);
                for (int j = 0; j < curcoords.size(); j++)
                    coords[3*ncn[i]*e+j] = curcoords[j];
            }
            element myelem(i,co);
            std::vector<std::vector<double>> tempsplit;
            myelem.fullsplit(mynumsplitrequested, tempsplit, coords);
            
            // Group with the existing splitcoords[p]:
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < tempsplit[j].size(); k++)
                    splitcoords[p][j][indexes[j]+k] = tempsplit[j][k];
                indexes[j] += tempsplit[j].size();
            }
        }
    }
    
    // All the info needed is ready. Clear and re-populate the nodes, elements and physical regions.
    mynodes = nodes();
    mydisjointregions = disjointregions();
    myphysicalregions = physicalregions(mydisjointregions);
    myelements = elements(mynodes, myphysicalregions, mydisjointregions);
    
    mynodes.setnumber(numnodes);
    std::vector<double>* nc = mynodes.getcoordinates();
    int nindex = 0;
    for (int p = 0; p < splitcoords.size(); p++)
    {
        physicalregion* curpr = myphysicalregions.get(prn[p]);
        
        for (int i = 0; i < 8; i++)
        {
            int ne = splitcoords[p][i].size()/3/ncn[i];
            for (int e = 0; e < ne; e++)
            {
                std::vector<int> nodelist = myalgorithm::getequallyspaced(nindex, 1, ncn[i]);
                int curelemindex = myelements.add(i, co, nodelist);
                curpr->addelement(i,curelemindex);
                // Add node coordinates:
                for (int k = 0; k < 3*ncn[i]; k++)
                    nc->at(3*nindex+k) = splitcoords[p][i][3*ncn[i]*e+k];
                    
                nindex += ncn[i];
            }
        }
    }
}

void rawmesh::readfromfile(std::string name)
{
    if (name == "gmshapi")
    {
        gmshinterface::readfromapi(mynodes, myelements, myphysicalregions);
        return;   
    }
    
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
    std::cout << "Use the GMSH API or the petsc mesh reader instead or use the GMSH .msh or Nastran .nas format." << std::endl;
    abort();
}

void rawmesh::writetofile(std::string name)
{
    if (name.length() >= 5 && name.compare(name.size()-4,4,".msh") == 0)
        gmshinterface::writetofile(name, mynodes, myelements, myphysicalregions, mydisjointregions);
    else
    {
        std::cout << "Error: file '" << name << "' has either no extension or it is not supported." << std::endl << "Currently supported: GMSH .msh" << std::endl;
        abort();
    }
}

void rawmesh::sortbybarycenters(int lasttypetoprocess)
{
    for (int elementtypenumber = 0; elementtypenumber <= lasttypetoprocess; elementtypenumber++)
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

void rawmesh::removeduplicates(int lasttypetoprocess)
{
    for (int elementtypenumber = 0; elementtypenumber <= lasttypetoprocess; elementtypenumber++)
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

void rawmesh::printcount(void)
{
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        if (myelements.count(elementtypenumber) == 0)
            continue;
            
        element elementobject(elementtypenumber);
        std::string elementname = elementobject.gettypenameconjugation(myelements.count(elementtypenumber));
        // For nodes there is no curvature order:
        if (elementtypenumber == 0)
            std::cout << "Extracted " << myelements.count(elementtypenumber) << " " << "nodes" << std::endl;
        else
            std::cout << "Extracted " << myelements.count(elementtypenumber) << " " << elementname << " with curvature order " << myelements.getcurvatureorder() << std::endl;
    }
}

rawmesh::rawmesh(void) : myelements(mynodes, myphysicalregions, mydisjointregions), myphysicalregions(mydisjointregions), myregiondefiner(mynodes, myelements, myphysicalregions) {}

nodes* rawmesh::getnodes(void) {return &mynodes;}
elements* rawmesh::getelements(void) {return &myelements;}
physicalregions* rawmesh::getphysicalregions(void) {return &myphysicalregions;}
disjointregions* rawmesh::getdisjointregions(void) {return &mydisjointregions;}
std::shared_ptr<ptracker> rawmesh::getptracker(void) {return myptracker;}
std::shared_ptr<htracker> rawmesh::gethtracker(void) {return myhtracker;}

std::shared_ptr<rawmesh> rawmesh::copy(void)
{
    std::shared_ptr<rawmesh> om(new rawmesh);

    *om = *this;
    
    om->mynodes = mynodes;
    om->mydisjointregions = mydisjointregions;
    myphysicalregions.copy(&(om->mydisjointregions), &(om->myphysicalregions));
    om->myelements = myelements.copy(&(om->mynodes), &(om->myphysicalregions), &(om->mydisjointregions));
    om->myregiondefiner = myregiondefiner.copy(&(om->mynodes), &(om->myelements), &(om->myphysicalregions));

    return om;
}

std::shared_ptr<rawmesh> rawmesh::getattarget(std::shared_ptr<ptracker> targetpt)
{
    // All calls from mathop::adapt immediately return (myptracker always equals targetpt).
    if (myptracker == targetpt)
        return shared_from_this();

    // Here we are in a rawfield syncing and the ptracker and htracker are not needed in 'om'.
    
    std::shared_ptr<rawmesh> om(new rawmesh);

    om->mynodes = mynodes;
    om->mydisjointregions = *(targetpt->getdisjointregions());
    om->myphysicalregions = physicalregions(om->mydisjointregions);
    om->myelements = myelements.copy(&(om->mynodes), &(om->myphysicalregions), &(om->mydisjointregions));

    (om->myelements).toptracker(myptracker, targetpt);
    
    return om;
}

void rawmesh::load(std::string name, int verbosity, bool legacyreader)
{
    // Do not call this when the mesh is already loaded!

    if (verbosity > 0)
    {
        if (name == "gmshapi")
            std::cout << "Loading mesh from GMSH API" << std::endl;
        else
            std::cout << "Loading mesh from file '" << name << "'" << std::endl;
    }
    
    wallclock loadtime;
    
    if (legacyreader || name == "gmshapi")
        readfromfile(name);
    else
    {
        petscmesh pmesh(name);
        pmesh.extract(mynodes, myelements, myphysicalregions);
    }
    
    splitmesh();
    mynodes.fixifaxisymmetric();
    
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
    
    mynumber = 0;
    
    myptracker = std::shared_ptr<ptracker>(new ptracker(myelements.count()));
    myptracker->updatedisjointregions(&mydisjointregions);
    
    myhtracker = std::shared_ptr<htracker>(new htracker(shared_from_this()));
    myhadaptedmesh = copy();
}

void rawmesh::load(bool mergeduplicates, std::vector<std::string> meshfiles, int verbosity)
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
    
    // Due to the shifting the bounds in regiondefiner will not be treated correctly:
    if (mergeduplicates == false && myregiondefiner.isanycoordinatedependentregiondefined())
    {
        std::cout << "Error in 'rawmesh' object: cannot define the requested region during loading when combining meshes without merging duplicates" << std::endl;
        abort();
    }
    
    load(flat, verbosity);
    
    // Shift back to the original position:
    if (mergeduplicates == false)
    {
        for (int i = 0; i < numfiles; i++)
        {
            shift(maxphysreg+1+i, -shiftvec[i],0,0);
            myhadaptedmesh->shift(maxphysreg+1+i, -shiftvec[i],0,0);
        }
    }
}

void rawmesh::load(std::vector<shape> inputshapes, int verbosity)
{
    // Do not call this when the mesh is already loaded!

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

    splitmesh();
    mynodes.fixifaxisymmetric();

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
    
    mynumber = 0;
    
    myptracker = std::shared_ptr<ptracker>(new ptracker(myelements.count()));
    myptracker->updatedisjointregions(&mydisjointregions);
    
    myhtracker = std::shared_ptr<htracker>(new htracker(shared_from_this()));
    myhadaptedmesh = copy();
}


void rawmesh::write(std::string name, int verbosity)
{
    if (verbosity > 0)
        std::cout << "Writing mesh to file '" << name << "'" << std::endl;
    
    wallclock writetime;
    
    writetofile(name);

    if (verbosity > 0)
        writetime.print("Time to write the mesh: ");
}

void rawmesh::split(int n)
{
    if (n < 0)
    {
        std::cout << "Error in 'mesh' object: number of splits cannot be negative" << std::endl;
        abort();
    }
    mynumsplitrequested += n;
}

std::vector<bool> rawmesh::isnodeinphysicalregion(int physreg)
{
    int numberofnodes = mynodes.count();
    
    std::vector<bool> output;

    if (physreg < 0)
        output = std::vector<bool>(numberofnodes, true);
    else
    {
        output = std::vector<bool>(numberofnodes, false);
    
        // Get only the disjoint regions with highest dimension elements:
        std::vector<int> selecteddisjregs = myphysicalregions.get(physreg)->getdisjointregions();
    
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            int disjreg = selecteddisjregs[i];
            int numelems = mydisjointregions.countelements(disjreg);
            int elemtypenum = mydisjointregions.getelementtypenumber(disjreg);
            int rangebegin = mydisjointregions.getrangebegin(disjreg);
            int curvatureorder = myelements.getcurvatureorder();
            
            element myelem(elemtypenum, curvatureorder);

            for (int e = 0; e < numelems; e++)
            {
                for (int n = 0; n < myelem.countcurvednodes(); n++)
                    output[myelements.getsubelement(0, elemtypenum, rangebegin+e, n)] = true;
            }
        }
    }
    
    return output;
}

void rawmesh::move(int physreg, expression u)
{
    int meshdim = getmeshdimension();
    if (u.countcolumns() != 1 || u.countrows() < meshdim || u.countrows() > 3)
    {
        std::cout << "Error in 'rawmesh' object: in 'move' expected a " << meshdim << "x1 expression" << std::endl;
        abort();
    }

    std::vector<double>* coords = mynodes.getcoordinates();
    
    u = u.resize(3,1);
    
    // To move each node only once:
    std::vector<bool> isnodemoved(mynodes.count(), false);
    
    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> selecteddisjregs;
    if (physreg >= 0)
        selecteddisjregs = (myphysicalregions.get(physreg))->getdisjointregions();
    else
        selecteddisjregs = mydisjointregions.getindim(meshdim);
        
    if (not(u.isharmonicone(selecteddisjregs)))
    {
        std::cout << "Error in 'rawmesh' object: the expression to move the mesh cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        abort();
    }
        
    // Send the disjoint regions with same element type numbers together:
    disjointregionselector mydisjregselector(selecteddisjregs, {});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
        int elementtypenumber = mydisjointregions.getelementtypenumber(mydisjregs[0]);

        // Get the reference coordinates corresponding to the nodes:
        lagrangeformfunction lff(elementtypenumber, myelements.getcurvatureorder(), {});
        std::vector<double> evaluationpoints = lff.getnodecoordinates();
        int ncn = evaluationpoints.size()/3;

        // Loop on all total orientations (if required):
        bool isorientationdependent = u.isvalueorientationdependent(mydisjregs);
        elementselector myselector(mydisjregs, isorientationdependent);
        do
        {            
            universe::allowreuse();

            // Compute the expression at the evaluation points:
            densematrix xval = u.getoperationinarray(0,0)->interpolate(myselector, evaluationpoints, NULL)[1][0];
            densematrix yval = u.getoperationinarray(1,0)->interpolate(myselector, evaluationpoints, NULL)[1][0];
            densematrix zval = u.getoperationinarray(2,0)->interpolate(myselector, evaluationpoints, NULL)[1][0];

            universe::forbidreuse();

            double* xvalptr = xval.getvalues();
            double* yvalptr = yval.getvalues();
            double* zvalptr = zval.getvalues();

            std::vector<int> elementnumbers = myselector.getelementnumbers();

            // Loop on all elements:
            for (int e = 0; e < elementnumbers.size(); e++)
            {
                for (int n = 0; n < ncn; n++)
                {
                    int curnode = myelements.getsubelement(0, elementtypenumber, elementnumbers[e], n);
                    if (isnodemoved[curnode] == false)
                    {
                        coords->at(3*curnode+0) += xvalptr[e*ncn+n];
                        coords->at(3*curnode+1) += yvalptr[e*ncn+n];
                        coords->at(3*curnode+2) += zvalptr[e*ncn+n];
                        
                        isnodemoved[curnode] = true;
                    }
                }
            }
        }
        while (myselector.next());
    }
    
    myelements.cleancoordinatedependentcontainers();
    
    mynodes.fixifaxisymmetric();
}

void rawmesh::shift(int physreg, double x, double y, double z)
{
    std::vector<double>* coords = mynodes.getcoordinates();

    std::vector<bool> isinsidereg = isnodeinphysicalregion(physreg);
    
    for (int n = 0; n < mynodes.count(); n++)
    {
        if (isinsidereg[n])
        {
            coords->at(3*n+0) += x;
            coords->at(3*n+1) += y;
            coords->at(3*n+2) += z;
        }
    }

    myelements.cleancoordinatedependentcontainers();
    
    mynodes.fixifaxisymmetric();
}

void rawmesh::rotate(int physreg, double ax, double ay, double az)
{
    std::vector<double>* coords = mynodes.getcoordinates();
    
    std::vector<bool> isinsidereg = isnodeinphysicalregion(physreg);
    
    std::vector<double> rotated = *coords;
    geotools::rotate(ax, ay, az, &rotated);
    
    for (int n = 0; n < mynodes.count(); n++)
    {
        if (isinsidereg[n])
        {
            coords->at(3*n+0) = rotated[3*n+0];
            coords->at(3*n+1) = rotated[3*n+1];
            coords->at(3*n+2) = rotated[3*n+2];
        }
    }
    
    myelements.cleancoordinatedependentcontainers();
    
    mynodes.fixifaxisymmetric();
}

void rawmesh::scale(int physreg, double x, double y, double z)
{
    std::vector<double>* coords = mynodes.getcoordinates();

    std::vector<bool> isinsidereg = isnodeinphysicalregion(physreg);
    
    for (int n = 0; n < mynodes.count(); n++)
    {
        if (isinsidereg[n])
        {
            coords->at(3*n+0) *= x;
            coords->at(3*n+1) *= y;
            coords->at(3*n+2) *= z;
        }
    }

    myelements.cleancoordinatedependentcontainers();
    
    mynodes.fixifaxisymmetric();
}

int rawmesh::getmeshdimension(void)
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

std::vector<int> rawmesh::getphysicalregionnumbers(int dim)
{
    if (dim < -1 || dim > 3)
    {
        std::cout << "Error in 'mesh' object: invalid input dimension '" << dim << "' in 'getphysicalregionnumbers'" << std::endl;
        abort();
    }

    return myphysicalregions.getallnumbers(dim);
}

void rawmesh::regionskin(int newphysreg, int physregtoskin)
{
    myregiondefiner.regionskin(newphysreg, physregtoskin);
}

void rawmesh::boxselection(int newphysreg, int physregtobox, int selecteddim, std::vector<double> boxlimit)
{
    myregiondefiner.boxselection(newphysreg, selecteddim, boxlimit, physregtobox);
}

void rawmesh::sphereselection(int newphysreg, int physregtosphere, int selecteddim, std::vector<double> centercoords, double radius)
{
    myregiondefiner.sphereselection(newphysreg, selecteddim, centercoords, radius, physregtosphere);
}

void rawmesh::regionexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude)
{
    myregiondefiner.regionexclusion(newphysreg, physregtoexcludefrom, physregstoexclude);
}


void rawmesh::add(std::shared_ptr<rawfield> inrawfield, expression criterion, std::vector<double> thresholds, std::vector<int> orders, double thresdown, double thresup, double mincritrange)
{
    int index = -1;
    for (int i = 0; i < mypadaptdata.size(); i++)
    {
        std::shared_ptr<rawfield> currawfield = (std::get<0>(mypadaptdata[i])).lock();
        if (currawfield.get() == inrawfield.get())
        {
            index = i;
            break;
        }
    }

    std::weak_ptr<rawfield> inweak = inrawfield;

    if (index != -1)
        mypadaptdata[index] = std::make_tuple(inweak, criterion, thresholds, orders, thresdown, thresup, mincritrange);
    else
        mypadaptdata.push_back(std::make_tuple(inweak, criterion, thresholds, orders, thresdown, thresup, mincritrange));
}

void rawmesh::remove(rawfield* inrawfield)
{
    // To delay criterion-field destruction to after 'mypadaptdata' has a valid state (after resize):
    std::vector<std::tuple<std::weak_ptr<rawfield>, expression,std::vector<double>,std::vector<int>,double,double,double>> pad = mypadaptdata;
    
    // Remove the pointed field and all expired fields:
    int curindex = 0;
    for (int i = 0; i < mypadaptdata.size(); i++)
    {
        std::weak_ptr<rawfield> curwp = std::get<0>(mypadaptdata[i]);
        if (not(curwp.expired()))
        {
            std::shared_ptr<rawfield> currawfield = curwp.lock();
            if (currawfield.get() != inrawfield)
            {
                if (curindex != i)
                    mypadaptdata[curindex] = mypadaptdata[i];
                curindex++;
            }
        }
    }
    mypadaptdata.resize(curindex);
}

bool rawmesh::adaptp(std::shared_ptr<rawmesh> critcalcrm, std::shared_ptr<ptracker> critcalcpt, bool washadapted, int verbosity)
{
    if (mypadaptdata.size() == 0)
        return false;

    double noisethreshold = 1e-8;
    int num = mypadaptdata.size();


    ///// Evaluate the criteria and the current interpolation orders:
    
    std::vector<vec> crits(num), oldords(num);
    
    int wholedomain = myphysicalregions.createunionofall();
    
    field one("one");
    for (int i = 0; i < num; i++)
    {
        std::shared_ptr<rawfield> curraw = (std::get<0>(mypadaptdata[i])).lock();
    
        expression critexpr = std::get<1>(mypadaptdata[i]);
        if (washadapted)
            critexpr = mathop::athp(critexpr, critcalcrm, critcalcpt);
        expression orderexpr = mathop::getfieldorder(field(curraw));
        if (washadapted)
            orderexpr = mathop::athp(orderexpr, critcalcrm, critcalcpt);
            
        formulation critaverage;
        critaverage += mathop::integral(wholedomain, -critexpr * mathop::tf(one) / mathop::getmeshsize(2) );
        critaverage.generaterhs();
        crits[i] = critaverage.rhs();
        
        formulation elorder;
        elorder += mathop::integral(wholedomain, -orderexpr * mathop::tf(one) / mathop::getmeshsize(2) );
        elorder.generaterhs();
        oldords[i] = elorder.rhs();
    }
    int totalnumelems = crits[0].size();

    // Move to densematrix container:
    std::vector<densematrix> critsmat(num), oldordsmat(num);
    std::vector<intdensematrix> newordsmat(num);
    intdensematrix ads(totalnumelems,1, 0, 1);
    for (int i = 0; i < num; i++)
    {
        critsmat[i] = crits[i].getvalues(ads);
        oldordsmat[i] = oldords[i].getvalues(ads);
        newordsmat[i] = intdensematrix(totalnumelems,1);
    }
    
    // Remove the above created region:
    myphysicalregions.remove({wholedomain}, false);

    // All vecs will have the same dofmanager:
    std::shared_ptr<dofmanager> dofmngr = crits[0].getpointer()->getdofmanager();
    dofmngr->selectfield(one.getpointer()->harmonic(1));
    crits = {}; oldords = {};


    ///// Calculate the new interpolation orders to use:
    
    int newmaxorder = 0;
    
    bool isidenticalorder = true;
    for (int i = 0; i < num; i++)
    {
        std::vector<double> thresholds = std::get<2>(mypadaptdata[i]);
        std::vector<int> orders = std::get<3>(mypadaptdata[i]);
        double thresdown = std::get<4>(mypadaptdata[i]) + noisethreshold;
        double thresup = std::get<5>(mypadaptdata[i]) + noisethreshold;
        double mincritrange = std::get<6>(mypadaptdata[i]);
        
        int minorder = *std::min_element(orders.begin(), orders.end());
    
        double* critsptr = critsmat[i].getvalues();
        double* oldordsptr = oldordsmat[i].getvalues();
        int* newordsptr = newordsmat[i].getvalues();
    
        std::vector<double> minmax = critsmat[i].minmax();
        double cmin = minmax[0];
        double cmax = minmax[1];
        
        double crange = cmax-cmin;
        
        // Convert the thresholds from % to criterion value:
        for (int th = 0; th < thresholds.size(); th++)
            thresholds[th] = cmin + thresholds[th] * crange;
        
        // In case the criterion range is not large enough the element orders are all set to minimum:
        if (crange < mincritrange)
        {
            for (int e = 0; e < totalnumelems; e++)
                newordsptr[e] = minorder;
            if (minorder > newmaxorder)
                newmaxorder = minorder;
        }
        else
        {
            for (int e = 0; e < totalnumelems; e++)
            {
                int oldord = (int)std::round(oldordsptr[e]);
                int interv = myalgorithm::findinterval(critsptr[e], thresholds);
                newordsptr[e] = orders[interv];
                
                // Check if the criterion is beyond the down/up change threshold.
                double intervsize = thresholds[interv+1]-thresholds[interv];
                // Bring back to upper interval?
                if (interv < thresholds.size()-2 && orders[interv+1] == oldord && critsptr[e] > thresholds[interv+1]-intervsize*thresdown)
                    newordsptr[e] = oldord;
                // Bring back to lower interval?
                if (interv > 0 && orders[interv-1] == oldord && critsptr[e] < thresholds[interv]+intervsize*thresup)
                    newordsptr[e] = oldord;
                
                if (newordsptr[e] > newmaxorder)
                    newmaxorder = newordsptr[e];
            }
        }
        
        for (int e = 0; e < totalnumelems; e++)
        {
            if (not(isidenticalorder) || std::abs(oldordsptr[e]-(double)newordsptr[e]) > 0.1)
            {
                isidenticalorder = false;
                break;
            }
        }
    }
    // Nothing to do if all new orders are identical to the old ones:
    if (isidenticalorder && not(washadapted))
    {
        if (verbosity > 0)
            std::cout << "Nothing to do for p-adaptation." << std::endl;
        return false;
    }
    
    
    ///// Add the elements to new physical regions:
    
    int lastphysregnum = myphysicalregions.getmaxphysicalregionnumber();
    int newlastpr = lastphysregnum;
    
    // newphysregsfororder[rawfield][order][elemtype].
    std::vector<std::vector<std::vector<physicalregion*>>> newphysregsfororder(num, std::vector<std::vector<physicalregion*>>(newmaxorder+1, std::vector<physicalregion*>(8, NULL)));
    
    for (int i = 0; i < num; i++)
    {
        int* newordsptr = newordsmat[i].getvalues();
        
        for (int d = 0; d < mydisjointregions.count(); d++)
        {
            if (dofmngr->isdefined(d, 0)) // There is only one shape fct!
            {
                int typenum = mydisjointregions.getelementtypenumber(d);
                int numelems = mydisjointregions.countelements(d);
                int rb = dofmngr->getrangebegin(d,0);
                int rbe = mydisjointregions.getrangebegin(d);
            
                for (int e = 0; e < numelems; e++)
                {
                    int curorder = newordsptr[rb+e];
        
                    if (newphysregsfororder[i][curorder][typenum] == NULL)
                    {
                        newphysregsfororder[i][curorder][typenum] = myphysicalregions.get(newlastpr+1);
                        newlastpr++;
                    }
                    
                    newphysregsfororder[i][curorder][typenum]->addelement(typenum, rbe+e);
                }
            }
        }
    }
    
    
    ///// Update the mesh:
    
    // The previous mesh tracker should not be touched:
    std::shared_ptr<ptracker> newptracker(new ptracker(myelements.count()));
    *newptracker = *myptracker;
    myptracker = newptracker;
    
    mydisjointregions.clear();
    
    myelements.definedisjointregions();

    std::vector<std::vector<int>> renumvec;
    myelements.reorderbydisjointregions(renumvec);
    
    myptracker->updaterenumbering(renumvec);
    myhtracker->renumbertransitions(renumvec);

    myelements.definedisjointregionsranges();
    
    // Define the physical regions based on the disjoint regions they contain:
    for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
        currentphysicalregion->definewithdisjointregions();
    }
    

    ///// New mesh version:
    
    myptracker->updatedisjointregions(&mydisjointregions);
    getoriginalmeshpointer()->mynumber++;
    mynumber = getoriginalmeshpointer()->mynumber;

    
    ///// Synchronize the rawfields and set their interpolation order:
    
    for (int i = 0; i < num; i++)
    {
        std::vector<int> curphysregsfororder = {};
        
        for (int o = 0; o < newphysregsfororder[i].size(); o++)
        {   
            for (int typenum = 0; typenum < 8; typenum++)
            {
                if (newphysregsfororder[i][o][typenum] != NULL)
                {
                    curphysregsfororder.push_back(newphysregsfororder[i][o][typenum]->getnumber());
                    curphysregsfororder.push_back(o);
                }
            }
        }
        
        std::shared_ptr<rawfield> curraw = (std::get<0>(mypadaptdata[i])).lock();
        
        // 'curphysregsfororder' MUST BE IN INCREASING ORDER.
        curraw->synchronize(curphysregsfororder);
    }
    

    ///// Remove the physical regions created:
    
    std::vector<int> prtoremove(newlastpr-lastphysregnum);
    for (int i = 0; i < prtoremove.size(); i++)
        prtoremove[i] = lastphysregnum+1+i;
    myphysicalregions.remove(prtoremove, true);
    
    myptracker->updatedisjointregions(&mydisjointregions);
    
    
    ///// Print p-adaptation summary:
    
    if (verbosity > 0)
        std::cout << "Adapted order of " << num << " field" << myalgorithm::getplurals(num) << "." << std::endl;
    if (verbosity > 1)
    {
        for (int i = 0; i < num; i++)
        {
            std::shared_ptr<rawfield> curraw = (std::get<0>(mypadaptdata[i])).lock();
            std::vector<int> numineachorder = newordsmat[i].countalloccurences(newmaxorder);

            curraw->print();

            for (int o = 0; o <= newmaxorder; o++)
                std::cout << " " << numineachorder[o];
            std::cout << std::endl;
        }
    }
    
    return true;
}

bool rawmesh::adapth(int verbosity)
{
    if (myhadaptdata.size() == 0)
        return false;
        
    wallclock clk;

    double noisethreshold = 1e-8;
    int meshdim = getmeshdimension();
    
    universe::mymesh = myhadaptedmesh;
        
    elements* elptr = universe::mymesh->getelements();
    disjointregions* drptr = universe::mymesh->getdisjointregions();
    physicalregions* prptr = universe::mymesh->getphysicalregions();
    
    std::shared_ptr<htracker> newhtracker(new htracker);
    *newhtracker = *(myhadaptedmesh->myhtracker);
    
    std::vector<int> leavesnumsplits;
    newhtracker->countsplits(leavesnumsplits);


    ///// Evaluate the criterion:

    int wholedomain = prptr->createunionofall();
    
    field one("one");
    
    formulation critaverage;
    critaverage += mathop::integral(wholedomain, -std::get<0>(myhadaptdata[0]) * mathop::tf(one) / mathop::getmeshsize(2) );
    critaverage.generaterhs();
    vec crit = critaverage.rhs();
    
    prptr->remove({wholedomain}, false);


    ///// Move to densematrix container:
    
    int totalnumelems = crit.size();

    intdensematrix ads(totalnumelems,1, 0, 1);
    densematrix critmat = crit.getvalues(ads);
    double* critptr = critmat.getvalues();

    std::shared_ptr<dofmanager> dofmngr = crit.getpointer()->getdofmanager();
    dofmngr->selectfield(one.getpointer()->harmonic(1));
    

    ///// Calculate the number of splits to use:
    
    std::vector<double> thresholds = std::get<1>(myhadaptdata[0]);
    std::vector<int> numsplits = std::get<2>(myhadaptdata[0]);
    double thresdown = std::get<3>(myhadaptdata[0]) + noisethreshold;
    double thresup = std::get<4>(myhadaptdata[0]) + noisethreshold;
    double mincritrange = std::get<5>(myhadaptdata[0]);
    
    int minnumsplits = *std::min_element(numsplits.begin(), numsplits.end());

    std::vector<double> minmax = critmat.minmax();
    double cmin = minmax[0];
    double cmax = minmax[1];
    
    double crange = cmax-cmin;
    
    // Convert the thresholds from % to criterion value:
    for (int th = 0; th < thresholds.size(); th++)
        thresholds[th] = cmin + thresholds[th] * crange;
    
    // Vector with +1 to split a leaf, 0 to keep as is and -1 to group:
    std::vector<int> vadapt(newhtracker->countleaves(), -1); // all grouped initially
    
    for (int d = 0; d < drptr->count(); d++)
    {
        if (dofmngr->isdefined(d, 0)) // There is only one shape fct!
        {
            int typenum = drptr->getelementtypenumber(d);
            int numelems = drptr->countelements(d);
            int rb = dofmngr->getrangebegin(d,0);
            int rbe = drptr->getrangebegin(d);
        
            for (int e = 0; e < numelems; e++)
            {
                int ln = newhtracker->getleafnumber(typenum, rbe+e);
                int oldnumsplits = leavesnumsplits[ln];
                
                int newnumsplits;
                    
                // In case the criterion range is not large enough the number of splits are all set to minimum:
                if (crange < mincritrange)
                    newnumsplits = minnumsplits;
                else
                {
                    double curcrit = critptr[rb+e];
                    
                    int interv = myalgorithm::findinterval(curcrit, thresholds);
                    newnumsplits = numsplits[interv];
                    
                    // Check if the criterion is beyond the down/up change threshold.
                    double intervsize = thresholds[interv+1]-thresholds[interv];
                    // Bring back to upper interval?
                    if (interv < thresholds.size()-2 && numsplits[interv+1] == oldnumsplits && curcrit > thresholds[interv+1]-intervsize*thresdown)
                        newnumsplits = oldnumsplits;
                    // Bring back to lower interval?
                    if (interv > 0 && numsplits[interv-1] == oldnumsplits && curcrit < thresholds[interv]+intervsize*thresup)
                        newnumsplits = oldnumsplits;
                }
                    
                // Split or keep unchanged. Multiple transition elements can share the same leaf!
                if (newnumsplits > oldnumsplits)
                    vadapt[ln] = 1;
                if (vadapt[ln] < 1 && newnumsplits == oldnumsplits)
                    vadapt[ln] = 0;
            }
        }
    }
             
        
    ///// Propagate the splits to guarantee at most a one delta between neighbouring elements:

    int maxnumsplits = newhtracker->getmaxdepth() + 1; // includes any new split request
    
    for (int ns = maxnumsplits; ns > 1; ns--)
    {
        // Update to how it will be actually treated:
        newhtracker->fix(vadapt);

        for (int d = 0; d < drptr->count(); d++)
        {
            if (dofmngr->isdefined(d, 0)) // There is only one shape fct!
            {
                int typenum = drptr->getelementtypenumber(d);
                int numelems = drptr->countelements(d);
                int rbe = drptr->getrangebegin(d);
                element curelem(typenum);
                int numedges = curelem.countedges();
            
                for (int e = 0; e < numelems; e++)
                {
                    int ln = newhtracker->getleafnumber(typenum, rbe+e);
                    int numsplits = leavesnumsplits[ln] + vadapt[ln];
                    
                    if (numsplits != ns)
                        continue;
                    
                    // Loop on every edge of the current element:
                    for (int en = 0; en < numedges; en++)
                    {
                        int currentedge = elptr->getsubelement(1, typenum, rbe+e, en);
                        // Get all cells touching the current edge:
                        std::vector<int> cellsonedge = elptr->getcellsonedge(currentedge);

                        for (int c = 0; c < cellsonedge.size()/2; c++)
                        {
                            int curcell = cellsonedge[2*c+1];
                            int celltype = cellsonedge[2*c+0];
                            
                            if (typenum == celltype && curcell == rbe+e)
                                continue;
                            
                            int curln = newhtracker->getleafnumber(celltype, curcell);
                            int neighbournumsplits = leavesnumsplits[curln] + vadapt[curln];

                            if (numsplits > neighbournumsplits+1)
                                vadapt[curln] += numsplits-neighbournumsplits-1;                     
                        }
                    }
                }
            }
        }
    }
        
        
    // Nothing to do if all new number of splits are identical to the old ones:
    newhtracker->fix(vadapt);
    bool isidentical = true;
    for (int i = 0; i < vadapt.size(); i++)
    {
        if (vadapt[i] != 0)
        {
            isidentical = false;
            break;
        }
    }
    if (isidentical)
    {
        if (verbosity > 0)
            std::cout << "Nothing to do for h-adaptation." << std::endl;
        return false;
    }
        
    newhtracker->adapt(vadapt);
    
    
    ///// Get the adapted element coordinates 'ac' and add to the mesh containers:
    
    myhadaptedmesh = std::shared_ptr<rawmesh>(new rawmesh);
    myhadaptedmesh->myhtracker = newhtracker;
    
    std::vector<std::vector<double>> ac;
    newhtracker->getadaptedcoordinates(ac, mynodes.getnoisethreshold());
    
    // Initialize nodes size:
    int numnodes = 0;
    for (int i = 0; i < 8; i++)
        numnodes += ac[i].size()/3;
    myhadaptedmesh->mynodes.setnumber(numnodes);
    std::vector<double>* nc = myhadaptedmesh->mynodes.getcoordinates();

    // Get the original element type and index for all leaves:
    std::vector<int> oet, oei;
    newhtracker->getoriginalelement(oet, oei);
    
    // Get the list of all physical regions in which each original element is:
    physicalregions* origpr = getphysicalregions();
    elements* origelems = getelements();

    std::vector<std::vector<int>> addresses(8, std::vector<int>(0));
    std::vector<std::vector<int>> prnums(8, std::vector<int>(0));
    for (int i = 0; i < 8; i++)
        origpr->inphysicalregions(i, origelems->count(i), addresses[i], prnums[i]);

    // Loop on all elements of max dimension in the adapted mesh:
    int ni = 0;
    int co = myelements.getcurvatureorder();
    for (int i = 0; i < 8; i++)
    {
        element myelem(i, co);
        int ncn = myelem.countcurvednodes();
        int numelems = ac[i].size()/ncn/3;
    
        for (int e = 0; e < numelems; e++)
        {
            int ln = newhtracker->getleafnumber(i, e);
            int origtype = oet[ln];
            int origindex = oei[ln];
            int numprs = addresses[origtype][origindex+1]-addresses[origtype][origindex];
            
            // Add nodes coordinates:
            for (int j = 0; j < 3*ncn; j++)
                nc->at(3*ni+j) = ac[i][3*ncn*e+j];
                
            // Add element:
            std::vector<int> nodes = myalgorithm::getequallyspaced(ni, 1, ncn);
            int elementindexincurrenttype = myhadaptedmesh->myelements.add(i, co, nodes);
            
            // Add the element to all required physical regions:
            for (int p = 0; p < numprs; p++)
                myhadaptedmesh->myphysicalregions.get(prnums[origtype][addresses[origtype][origindex]+p])->addelement(i, elementindexincurrenttype);
            
            ni += ncn;
        }
    }
    
    myhadaptedmesh->mynodes.fixifaxisymmetric();


    ///// Process the mesh:
        
    myhadaptedmesh->myelements.explode();

    std::vector<int> lasttype = {-1,0,1,3};
    myhadaptedmesh->sortbybarycenters(lasttype[meshdim]);
    myhadaptedmesh->removeduplicates(lasttype[meshdim]);
    
    
    ///// Add the lower dimension physical regions:
    
    std::vector<int> dims = {0,1,2,2,3,3,3,3};
    // Check if there is any node/edge/face physical region for speedup:
    std::vector<bool> isanypratdim(4, false);
    for (int i = 0; i < 8; i++)
        isanypratdim[dims[i]] = (isanypratdim[dims[i]] || prnums[i].size() > 0);
    
    // Loop on all transition elements:
    for (int i = 0; i < 8; i++)
    {
        if (dims[i] != meshdim)
            continue;
        
        std::vector<int> on, oe, of;
        for (int j = 0; j < myhadaptedmesh->myelements.count(i); j++)
        {
            int origtype, origelemnum;
            newhtracker->atoriginal(i, j, origtype, origelemnum, on, oe, of);
            
            if (isanypratdim[0])
            {
                for (int n = 0; n < on.size(); n++)
                {
                    if (on[n] != -1)
                    {
                        int curnum = myhadaptedmesh->myelements.getsubelement(0, i, j, n);
                        int curorignum = origelems->getsubelement(0, origtype, origelemnum, on[n]);
                        // Loop on all physical regions for this subelement:
                        int numpr = addresses[0][curorignum+1]-addresses[0][curorignum];
                        for (int p = 0; p < numpr; p++)
                            myhadaptedmesh->myphysicalregions.get( prnums[0][addresses[0][curorignum]+p] )->addelement(0,curnum);
                    }
                }
            }
            if (isanypratdim[1])
            {
                for (int e = 0; e < oe.size(); e++)
                {
                    if (oe[e] != -1)
                    {
                        int curnum = myhadaptedmesh->myelements.getsubelement(1, i, j, e);
                        int curorignum = origelems->getsubelement(1, origtype, origelemnum, oe[e]);
                        // Loop on all physical regions for this subelement:
                        int numpr = addresses[1][curorignum+1]-addresses[1][curorignum];
                        for (int p = 0; p < numpr; p++)
                            myhadaptedmesh->myphysicalregions.get( prnums[1][addresses[1][curorignum]+p] )->addelement(1,curnum);
                    }
                }
            }
            if (isanypratdim[2])
            {
                for (int f = 0; f < of.size(); f++)
                {
                    if (of[f] != -1)
                    {
                        // FOR HEXAHEDRA, PRISMS AND PYRAMIDS THIS IS INCORRECT SINCE THEIR FACES ARE NOT ONLY TRIANGULAR
                        // DO NOT FORGET THAT THE TRANSITION AND ORIGINAL ELEMENTS MIGHT BE OF DIFFERENT TYPES AND THE FACES AS WELL!
                        int curnum = myhadaptedmesh->myelements.getsubelement(2, i, j, f);
                        int curorignum = origelems->getsubelement(2, origtype, origelemnum, of[f]);
                        // Loop on all physical regions for this subelement:
                        int numpr = addresses[2][curorignum+1]-addresses[2][curorignum];
                        for (int p = 0; p < numpr; p++)
                            myhadaptedmesh->myphysicalregions.get( prnums[2][addresses[2][curorignum]+p] )->addelement(2,curnum);
                    }
                }
            }
        }
    }
    
    // Remove the duplicated elements in every physical region:
    for (int physregindex = 0; physregindex < myhadaptedmesh->myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myhadaptedmesh->myphysicalregions.getatindex(physregindex);
        currentphysicalregion->removeduplicatedelements();
    }
    
    
    ///// Continue processing the mesh:
    
    myhadaptedmesh->myelements.definedisjointregions();
    
    std::vector<std::vector<int>> renumbydr;
    myhadaptedmesh->myelements.reorderbydisjointregions(renumbydr);
    newhtracker->renumbertransitions(renumbydr);
    
    myhadaptedmesh->myelements.definedisjointregionsranges();
    
    // Define the physical regions based on the disjoint regions they contain:
    for (int physregindex = 0; physregindex < myhadaptedmesh->myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myhadaptedmesh->myphysicalregions.getatindex(physregindex);
        currentphysicalregion->definewithdisjointregions();
    }
    
    myhadaptedmesh->myelements.orient();

    
    // Optional output:
    if (verbosity > 0)
    {
        std::string out = "";
        for (int i = 0; i < 8; i++)
        {
            element myelem(i);
            int curnum = myhadaptedmesh->myelements.count(i);
            if (dims[i] == meshdim && curnum > 0)
                out = out + std::to_string(curnum) + " " + myelem.gettypenameconjugation(curnum) + " + ";
        }
        out.resize(out.size()-3);
    
        clk.print("Adapted to " + out + " (" + std::to_string(newhtracker->getmaxdepth()) + ") in");
    }
    
    
    ///// New mesh version:
    
    mynumber++;
    myhadaptedmesh->mynumber = mynumber;
    
    myhadaptedmesh->myptracker = std::shared_ptr<ptracker>(new ptracker(myhadaptedmesh->myelements.count()));
    myhadaptedmesh->myptracker->updatedisjointregions(&(myhadaptedmesh->mydisjointregions));
    myhadaptedmesh->mypadaptdata = universe::mymesh->mypadaptdata;
    universe::mymesh->mypadaptdata = {};
    
    
    ///// Send mesh to universe:
    universe::mymesh = myhadaptedmesh;
    
    
    return true;
}

void rawmesh::setadaptivity(expression criterion, int lownumsplits, int highnumsplits, double thresdown, double thresup, double mincritrange)
{
    if (lownumsplits < 0)
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' cannot use negative minimum number of splits " << lownumsplits << std::endl;
        abort();   
    }
    if (highnumsplits < lownumsplits)
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' the minimum number of splits cannot be larger than the maximum" << std::endl;
        abort();   
    }

    int cnt = highnumsplits-lownumsplits+1;
    std::vector<double> thresholds(cnt+1);
    thresholds[0] = 0;
    for (int i = 1; i < thresholds.size()-1; i++)
        thresholds[i] = 1.0/cnt * i;
    thresholds[thresholds.size()-1] = 1.0;

    std::vector<int> numsplits(cnt);
    for (int i = 0; i < cnt; i++)
        numsplits[i] = lownumsplits + i;
        
    setadaptivity(criterion, thresholds, numsplits, thresdown, thresup, mincritrange);
}

void rawmesh::setadaptivity(expression criterion, std::vector<double> thresholds, std::vector<int> numsplits, double thresdown, double thresup, double mincritrange)
{
    double noiselevel = 1e-8;

    if (not(criterion.isscalar()))
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' expected a scalar criterion for mesh adaptivity" << std::endl;
        abort();   
    }
    if (numsplits.size() == 0)
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' number of splits vector cannot be empty" << std::endl;
        abort();   
    }
    if (thresholds.size() != numsplits.size()+1)
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' expected a threshold vector of size " << numsplits.size()+1 << " for the number of splits vector of size " << numsplits.size() << " provided (thresholds 0.0 and 1.0 must be included)" << std::endl;
        abort();   
    }
    for (int i = 0; i < thresholds.size(); i++)
    {
        if (thresholds[i] < 0.0-noiselevel || thresholds[i] > 1.0+noiselevel)
        {
            std::cout << "Error in 'rawmesh' object: in 'setadaptivity' thresholds must be between 0.0 and 1.0" << std::endl;
            abort();   
        }
    }
    for (int i = 1; i < thresholds.size(); i++)
    {
        if (thresholds[i] < thresholds[i-1]+noiselevel)
        {
            std::cout << "Error in 'rawmesh' object: in 'setadaptivity' expecting increasing thresholds" << std::endl;
            abort();   
        }
    }
    for (int i = 1; i < numsplits.size(); i++)
    {
        if (numsplits[i] <= numsplits[i-1])
        {
            std::cout << "Error in 'rawmesh' object: in 'setadaptivity' expecting increasing number of splits" << std::endl;
            abort();   
        }
    }
    if (thresholds[0] > noiselevel || thresholds[thresholds.size()-1] < 1.0-noiselevel)
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' thresholds 0.0 and 1.0 must be included" << std::endl;
        abort();   
    }
    for (int i = 0; i < numsplits.size(); i++)
    {
        if (numsplits[i] < 0)
        {
            std::cout << "Error in 'rawmesh' object: in 'setadaptivity' number of splits cannot be negative" << std::endl;
            abort();   
        }
    }
    if (thresdown < 0.0-noiselevel || thresdown > 1.0+noiselevel || thresup < 0.0-noiselevel || thresup > 1.0+noiselevel)
    {
        std::cout << "Error in 'rawmesh' object: in 'setadaptivity' thresholds must be between 0.0 and 1.0" << std::endl;
        abort();   
    }
    
    // The criterion cannot be multiharmonic:
    std::vector<int> alldisjregs(universe::mymesh->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(criterion.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'rawmesh' object: cannot have a multiharmonic criterion for h-adaptivity" << std::endl;
        abort();
    }
    
    myhadaptdata = {std::make_tuple(criterion, thresholds, numsplits, thresdown, thresup, mincritrange)};   
}

void rawmesh::writewithdisjointregions(std::string name)
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

void rawmesh::printphysicalregions(void)
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

void rawmesh::printdisjointregions(void)
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

void rawmesh::printelementsinphysicalregions(bool isdebug)
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

std::shared_ptr<rawmesh> rawmesh::gethadaptedpointer(void)
{
    if (myhadaptedmesh == NULL)
        return shared_from_this();
    else
        return myhadaptedmesh;
}

std::shared_ptr<rawmesh> rawmesh::getoriginalmeshpointer(void)
{
    return myhtracker->getoriginalmesh();
}

