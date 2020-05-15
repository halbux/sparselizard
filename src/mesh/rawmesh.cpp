#include "rawmesh.h"


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

void rawmesh::sortbybarycenters(void)
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

void rawmesh::removeduplicates(void)
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

void rawmesh::load(std::string name, int verbosity, bool legacyreader)
{
    ///// Reset all memory of the rawmesh object:
    mynodes = nodes();
    mydisjointregions = disjointregions();
    myphysicalregions = physicalregions(mydisjointregions);
    myelements = elements(mynodes, myphysicalregions, mydisjointregions);
    ///// Memory is reset


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
    myptracker = std::shared_ptr<ptracker>(new ptracker);
    myptracker->updatedisjointregions(&mydisjointregions);
    mypadaptdata = {};
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
    
    this->load(flat, verbosity);
    
    // Shift back to the original position:
    if (mergeduplicates == false)
    {
        for (int i = 0; i < numfiles; i++)
            this->shift(maxphysreg+1+i, -shiftvec[i],0,0);
    }
    
    mynumber = 0;
    myptracker = std::shared_ptr<ptracker>(new ptracker);
    myptracker->updatedisjointregions(&mydisjointregions);
    mypadaptdata = {};
}

void rawmesh::load(std::vector<shape> inputshapes, int verbosity)
{
    ///// Reset all memory of the rawmesh object:
    mynodes = nodes();
    mydisjointregions = disjointregions();
    myphysicalregions = physicalregions(mydisjointregions);
    myelements = elements(mynodes, myphysicalregions, mydisjointregions);
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

    splitmesh();

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
    myptracker = std::shared_ptr<ptracker>(new ptracker);
    myptracker->updatedisjointregions(&mydisjointregions);
    mypadaptdata = {};
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

void rawmesh::shift(int physreg, double x, double y, double z)
{
    mynodes.shift(physreg, x, y, z);
    myelements.shift(physreg, x, y, z);
}

void rawmesh::shift(double x, double y, double z)
{
    mynodes.shift(-1, x, y, z);
    myelements.shift(-1, x, y, z);
}

void rawmesh::rotate(int physreg, double ax, double ay, double az)
{
    mynodes.rotate(physreg, ax, ay, az);
    myelements.rotate(physreg, ax, ay, az);
}

void rawmesh::rotate(double ax, double ay, double az)
{
    mynodes.rotate(-1, ax, ay, az);
    myelements.rotate(-1, ax, ay, az);
}

void rawmesh::scale(int physreg, double x, double y, double z)
{
    mynodes.scale(physreg, x, y, z);
    myelements.scale(physreg, x, y, z);
}

void rawmesh::scale(double x, double y, double z)
{
    mynodes.scale(-1, x, y, z);
    myelements.scale(-1, x, y, z);
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

void rawmesh::adaptp(void)
{
    if (universe::ispadaptallowed == false)
        return;

    double noisethreshold = 1e-8;
    int num = mypadaptdata.size();


    ///// Evaluate the criteria and the current interpolation orders:
    
    std::vector<vec> crits(num), oldords(num);
    
    int wholedomain = mathop::regionall();
    
    field one("one");
    for (int i = 0; i < num; i++)
    {
        std::shared_ptr<rawfield> curraw = (std::get<0>(mypadaptdata[i])).lock();
    
        formulation critaverage;
        critaverage += mathop::integral(wholedomain, -std::get<1>(mypadaptdata[i]) * mathop::tf(one) / mathop::getmeshsize(2) );
        critaverage.generaterhs();
        crits[i] = critaverage.rhs();
        
        formulation elorder;
        elorder += mathop::integral(wholedomain, -mathop::getfieldorder(field(curraw)) * mathop::tf(one) / mathop::getmeshsize(2) );
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
    if (isidenticalorder)
        return;
        
    
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
    std::shared_ptr<ptracker> newptracker(new ptracker);
    *newptracker = *myptracker;
    myptracker = newptracker;
    
    mydisjointregions.clear();
    
    myelements.definedisjointregions();

    std::vector<std::vector<int>> renumvec;
    myelements.reorderbydisjointregions(renumvec);
    
    myptracker->updaterenumbering(renumvec);

    myelements.definedisjointregionsranges();
    
    // Define the physical regions based on the disjoint regions they contain:
    for (int physregindex = 0; physregindex < myphysicalregions.count(); physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions.getatindex(physregindex);
        currentphysicalregion->definewithdisjointregions();
    }

    
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
    
    
    ///// New mesh version:
    
    myptracker->updatedisjointregions(&mydisjointregions);
    mynumber++;
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

std::shared_ptr<rawmesh> rawmesh::getpointer(void)
{
    return shared_from_this();
}

