#include "elements.h"
#include "geotools.h"


elements::elements(nodes& inputnodes, physicalregions& inputphysicalregions, disjointregions& inputdisjointregions)
{
    mynodes = &inputnodes;
    myphysicalregions = &inputphysicalregions;
    mydisjointregions = &inputdisjointregions;
}

nodes* elements::getnodes(void) {return mynodes;}
physicalregions* elements::getphysicalregions(void) {return myphysicalregions;}
disjointregions* elements::getdisjointregions(void) {return mydisjointregions;}

int elements::add(int elementtypenumber, int curvatureorder, std::vector<int>& nodelist)
{
    // Point elements have a number equal to the node defining them
    // and should not be added to 'subelementsinelements'.
    if (elementtypenumber == 0)
        return nodelist[0];

    // Define the curvature order for the first time:
    if (mycurvatureorder == -1)
    {
        mycurvatureorder = curvatureorder;
        // Define as well the following quantities now that the curvature order is known:
        for (int i = 0; i <= 7; i++)
        {
            element myelement(i,curvatureorder);
            
            numberofsubelementsineveryelement[i][0] = myelement.countcurvednodes();
            numberofsubelementsineveryelement[i][1] = myelement.countedges();
            numberofsubelementsineveryelement[i][2] = myelement.counttriangularfaces();
            numberofsubelementsineveryelement[i][3] = myelement.countquadrangularfaces();
        }
    }
    if (mycurvatureorder != curvatureorder)
    {
        std::cout << "Error in 'elements' object: the mesh can contain only a single curvature order." << std::endl;
        abort();
    }
        
    subelementsinelements[elementtypenumber][0].insert(subelementsinelements[elementtypenumber][0].end(), nodelist.begin(), nodelist.end());
    return subelementsinelements[elementtypenumber][0].size()/nodelist.size() - 1;
}

void elements::cleancoordinatedependentcontainers(void)
{
    barycenters = std::vector<std::vector<double>>(8, std::vector<double>(0));
    sphereradius = std::vector<std::vector<double>>(8, std::vector<double>(0));
    boxdimensions = std::vector<std::vector<double>>(8, std::vector<double>(0));
}

int elements::getsubelement(int subelementtypenumber, int elementtypenumber, int elementnumber, int subelementindex)
{
    if (elementtypenumber == subelementtypenumber)
        return elementnumber;
    else
        return subelementsinelements[elementtypenumber][subelementtypenumber][elementnumber*numberofsubelementsineveryelement[elementtypenumber][subelementtypenumber]+subelementindex];
}

int elements::getdisjointregion(int elementtypenumber, int elementnumber, bool errorifnegative)
{
    int output = indisjointregion[elementtypenumber][elementnumber];
    if (output >= 0 || not(errorifnegative))
        return output;
    else
    {
        std::cout << "Error in 'elements' object: returned a negative disjoint region number" << std::endl;
        std::cout << "Possible reason: too bad quality elements obtained when h-adapting a curved mesh (some identical edges might not have been merged and therefore AMR has failed)" << std::endl;
        abort();
    }
}

int elements::gettotalorientation(int elementtypenumber, int elementnumber)
{
    if (elementtypenumber != 0)
        return totalorientations[elementtypenumber][elementnumber];
    else
        return 0;
}


int elements::count(int elementtypenumber)
{
    if (elementtypenumber == 0)
        return mynodes->count();
    
    if (subelementsinelements[elementtypenumber][0].size() == 0)
        return 0;
    
    return subelementsinelements[elementtypenumber][0].size()/numberofsubelementsineveryelement[elementtypenumber][0];
}

int elements::countindim(int dim)
{
    int output = 0;
    
    std::vector<int> indim = {0,1,2,2,3,3,3,3};
    for (int i = 0; i < 8; i++)
    {
        if (indim[i] == dim)
            output += count(i);
    }
    return output;
}

std::vector<int> elements::count(void)
{
    std::vector<int> output(8);
    for (int i = 0; i < 8; i++)
        output[i] = count(i);

    return output;
}

int elements::getcurvatureorder(void)
{
    return mycurvatureorder;
}

void elements::populateedgesatnodes(void)
{
    // Get the number of nodes:
    int numnodes = count(0);
    // Get the number of edges:
    int numedges = count(1);
    // Get the number of curved nodes per edge:
    int numcurvednodes = subelementsinelements[1][0].size()/numedges;


    // Preallocate vectors:
    adressedgesatnodes = std::vector<int>(numnodes,0);
    edgesatnodes = std::vector<int>(numedges*2);


    // Vector to count the number of edges touching every node:
    std::vector<int> numedgesonnode(numnodes,0);

    // Loop on all edges:
    for (int e = 0; e < numedges; e++)
    {
        // Corner nodes in current edge:
        int edgefirstnode = subelementsinelements[1][0][e*numcurvednodes+0];
        int edgelastnode = subelementsinelements[1][0][e*numcurvednodes+1];

        numedgesonnode[edgefirstnode]++;
        numedgesonnode[edgelastnode]++;
    }

    for (int n = 1; n < numnodes; n++)
        adressedgesatnodes[n] = adressedgesatnodes[n-1] + numedgesonnode[n-1];

    std::vector<int> currentedgeinnode(numnodes,0);    
    for (int e = 0; e < numedges; e++)
    {
        // Corner nodes in current edge:
        int edgefirstnode = subelementsinelements[1][0][e*numcurvednodes+0];
        int edgelastnode = subelementsinelements[1][0][e*numcurvednodes+1];

        edgesatnodes[adressedgesatnodes[edgefirstnode]+currentedgeinnode[edgefirstnode]] = e;
        currentedgeinnode[edgefirstnode]++;
        edgesatnodes[adressedgesatnodes[edgelastnode]+currentedgeinnode[edgelastnode]] = e;
        currentedgeinnode[edgelastnode]++;
    }
}

int elements::countedgesonnode(int nodenumber)
{
    if (edgesatnodes.size() == 0)
        populateedgesatnodes();

    if (nodenumber+1 < adressedgesatnodes.size())
        return adressedgesatnodes[nodenumber+1]-adressedgesatnodes[nodenumber];
    else
        return edgesatnodes.size()-adressedgesatnodes[nodenumber];
}

std::vector<int> elements::getedgesonnode(int nodenumber)
{
    int numedgesatnode = countedgesonnode(nodenumber);

    std::vector<int> output(numedgesatnode);
    for (int i = 0; i < numedgesatnode; i++)
        output[i] = edgesatnodes[adressedgesatnodes[nodenumber]+i];

    return output;
}

void elements::populatecellsatedges(void)
{
    int celldim = getdimension();

    // Get the number of edges:
    int numedges = count(1);
    // Get the number of cells of each type and the number of edges in each cell:
    std::vector<int> numcells(8,0);
    std::vector<int> ne(8);
    int prealloc = 0;
    for (int i = 0; i < 8; i++)
    {
        element myelem(i);
        if (myelem.getelementdimension() == celldim)
        {
            numcells[i] = count(i);
            ne[i] = myelem.countedges();
            prealloc += numcells[i]*ne[i];   
        }
    }


    // Preallocate vectors:
    adresscellsatedges = std::vector<int>(numedges,0);
    cellsatedges = std::vector<int>(2*prealloc);


    // Vector to count the number of cells touching every edge:
    std::vector<int> numcellsonedge(numedges,0);

    // Loop on all cells:
    for (int i = 0; i < 8; i++)
    {
        for (int c = 0; c < numcells[i]; c++)
        {
            for (int e = 0; e < ne[i]; e++)
                numcellsonedge[getsubelement(1,i,c,e)]++;
        }
    }

    for (int e = 1; e < numedges; e++)
        adresscellsatedges[e] = adresscellsatedges[e-1] + 2*numcellsonedge[e-1];

    std::vector<int> currentcellinedge(numedges,0); 
    for (int i = 0; i < 8; i++)
    {
        for (int c = 0; c < numcells[i]; c++)
        {
            for (int e = 0; e < ne[i]; e++)
            {
                int curedge = getsubelement(1,i,c,e);
                cellsatedges[adresscellsatedges[curedge]+2*currentcellinedge[curedge]+0] = i;
                cellsatedges[adresscellsatedges[curedge]+2*currentcellinedge[curedge]+1] = c;
                
                currentcellinedge[curedge]++;
            }
        }
    }
}

int elements::countcellsonedge(int edgenumber)
{
    if (cellsatedges.size() == 0)
        populatecellsatedges();

    if (edgenumber+1 < adresscellsatedges.size())
        return (adresscellsatedges[edgenumber+1]-adresscellsatedges[edgenumber])/2;
    else
        return (cellsatedges.size()-adresscellsatedges[edgenumber])/2;
}

std::vector<int> elements::getcellsonedge(int edgenumber)
{
    int numcellsatedge = countcellsonedge(edgenumber);

    std::vector<int> output(2*numcellsatedge);
    for (int i = 0; i < 2*numcellsatedge; i++)
        output[i] = cellsatedges[adresscellsatedges[edgenumber]+i];

    return output;
}

std::vector<double> elements::getnodecoordinates(int elementtypenumber, int elementnumber, int xyz)
{
    std::vector<double>* nodecoordinates = mynodes->getcoordinates();
    
    if (elementtypenumber == 0)
        return {nodecoordinates->at(3*elementnumber+xyz)};
    
    element myelement(elementtypenumber, mycurvatureorder);
    int curvednumberofnodes = myelement.countcurvednodes();
    
    std::vector<double> nodecoords(curvednumberofnodes);

    for (int node = 0; node < curvednumberofnodes; node++)
        nodecoords[node] = nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elementnumber*curvednumberofnodes+node]+xyz);

    return nodecoords;
}

std::vector<double> elements::getnodecoordinates(int elementtypenumber, int elementnumber)
{
    std::vector<double>* nodecoordinates = mynodes->getcoordinates();
    
    if (elementtypenumber == 0)
        return {nodecoordinates->at(3*elementnumber+0), nodecoordinates->at(3*elementnumber+1), nodecoordinates->at(3*elementnumber+2)};
    
    element myelement(elementtypenumber, mycurvatureorder);
    int curvednumberofnodes = myelement.countcurvednodes();
    
    std::vector<double> nodecoords(3*curvednumberofnodes);

    for (int node = 0; node < curvednumberofnodes; node++)
    {
        int curnode = subelementsinelements[elementtypenumber][0][elementnumber*curvednumberofnodes+node];
        nodecoords[3*node+0] = nodecoordinates->at(3*curnode+0);
        nodecoords[3*node+1] = nodecoordinates->at(3*curnode+1);
        nodecoords[3*node+2] = nodecoordinates->at(3*curnode+2);
    }

    return nodecoords;
}

void elements::getrefcoordsondisjregs(int origintype, std::vector<int>& elems, std::vector<double>& refcoords, std::vector<int> targetdisjregs, std::vector<int>& targetelems, std::vector<double>& targetrefcoords)
{
    std::vector<int> renumtoelems(count(origintype), -1);
    std::vector<bool> isprocessed(elems.size(), false);
    for (int i = 0; i < elems.size(); i++)
        renumtoelems[elems[i]] = i;

    int numrc = refcoords.size()/3;
    int numpoints = elems.size() * numrc;
    element mysubelement(origintype);
    int subnumnodes = mysubelement.countnodes();
    std::vector<double> cornerrcs(3*subnumnodes);
    
    targetelems = std::vector<int>(2*numpoints, -1);
    targetrefcoords = std::vector<double>(3*numpoints);
    
    std::vector<int> nodeindexincurelem(count(0));

    for (int i = 0; i < targetdisjregs.size(); i++)
    {
        int d = targetdisjregs[i];
        int tn = mydisjointregions->getelementtypenumber(d);
        int ne = mydisjointregions->countelements(d);
        int rb = mydisjointregions->getrangebegin(d);
        
        element myelement(tn);
        int numnodes = myelement.countnodes();
        int numsubelems = myelement.counttype(origintype);
        
        lagrangeformfunction lff(tn, 1, {});
        std::vector<double> rcs = lff.getnodecoordinates();
        
        // Loop on all candidate target elements:
        for (int e = 0; e < ne; e++)
        {
            for (int n = 0; n < numnodes; n++)
                nodeindexincurelem[getsubelement(0, tn, rb+e, n)] = n;
        
            // Loop on all 'origintype' subelements in the target:
            for (int s = 0; s < numsubelems; s++)
            {
                int cursubelem = getsubelement(origintype, tn, rb+e, s);
                int origindex = renumtoelems[cursubelem];
                
                if (origindex != -1 && isprocessed[origindex] == false)
                {
                    isprocessed[origindex] = true;
        
                    for (int n = 0; n < subnumnodes; n++)
                    {
                        int curnode = getsubelement(0, origintype, cursubelem, n);

                        cornerrcs[3*n+0] = rcs[3*nodeindexincurelem[curnode]+0];
                        cornerrcs[3*n+1] = rcs[3*nodeindexincurelem[curnode]+1];
                        cornerrcs[3*n+2] = rcs[3*nodeindexincurelem[curnode]+2];
                    }
                        
                    std::vector<double> rcsintarget = mysubelement.calculatecoordinates(refcoords, cornerrcs, 0, (origintype == 0));
                    
                    for (int r = 0; r < numrc; r++)
                    {
                        targetelems[2*origindex*numrc+2*r+0] = tn;
                        targetelems[2*origindex*numrc+2*r+1] = rb+e;
                        
                        targetrefcoords[3*origindex*numrc+3*r+0] = rcsintarget[3*r+0];
                        targetrefcoords[3*origindex*numrc+3*r+1] = rcsintarget[3*r+1];
                        targetrefcoords[3*origindex*numrc+3*r+2] = rcsintarget[3*r+2];
                    }
                }
            }
        }
    }
}

std::vector<double>* elements::getbarycenters(int elementtypenumber)
{
    // If not yet populated for the element type:
    if (barycenters[elementtypenumber].size() == 0)
        barycenters[elementtypenumber] = computebarycenters(elementtypenumber);
    
    return &(barycenters[elementtypenumber]);
}

std::vector<double>* elements::getsphereradius(int elementtypenumber)
{
    std::vector<double>* mybarys = getbarycenters(elementtypenumber);

    // If not yet populated for the element type:
    if (sphereradius[elementtypenumber].size() == 0)
    {
        sphereradius[elementtypenumber].resize(count(elementtypenumber));
    
        for (int i = 0; i < count(elementtypenumber); i++)
        {
            double maxdist = 0;
        
            std::vector<double> xnodes = getnodecoordinates(elementtypenumber, i, 0);
            std::vector<double> ynodes = getnodecoordinates(elementtypenumber, i, 1);
            std::vector<double> znodes = getnodecoordinates(elementtypenumber, i, 2);
            
            for (int j = 0; j < xnodes.size(); j++)
            {
                double curdist = std::sqrt( std::pow(mybarys->at(3*i+0)-xnodes[j], 2) + std::pow(mybarys->at(3*i+1)-ynodes[j], 2) + std::pow(mybarys->at(3*i+2)-znodes[j], 2) );
                if (curdist > maxdist)
                    maxdist = curdist;
            }
            
            sphereradius[elementtypenumber][i] = maxdist;
        }
    }
    
    return &(sphereradius[elementtypenumber]);
}

std::vector<double>* elements::getboxdimensions(int elementtypenumber)
{
    std::vector<double>* mybarys = getbarycenters(elementtypenumber);
    
    // If not yet populated for the element type:
    if (boxdimensions[elementtypenumber].size() == 0)
    {
        boxdimensions[elementtypenumber].resize(3*count(elementtypenumber));
    
        for (int i = 0; i < count(elementtypenumber); i++)
        {
            for (int c = 0; c < 3; c++)
            {
                std::vector<double> curnodescoords = getnodecoordinates(elementtypenumber, i, c);
                
                double maxdist = 0;
                for (int j = 0; j < curnodescoords.size(); j++)
                {
                    double curdist = std::abs(mybarys->at(3*i+c)-curnodescoords[j]);
                    if (curdist > maxdist)
                        maxdist = curdist;
                }
                boxdimensions[elementtypenumber][3*i+c] = maxdist;
            }
        }
    }
    
    return &(boxdimensions[elementtypenumber]);
}

std::vector<double> elements::getnormal(int elementtypenumber, int elementnumber)
{
    std::vector<double>* nodecoordinates = mynodes->getcoordinates();
    
    element myelement(elementtypenumber, mycurvatureorder);
    int curvednumberofnodes = myelement.countcurvednodes();
    
    int node0 = subelementsinelements[elementtypenumber][0][elementnumber*curvednumberofnodes+0];
    int node1 = subelementsinelements[elementtypenumber][0][elementnumber*curvednumberofnodes+1];
    int node2 = subelementsinelements[elementtypenumber][0][elementnumber*curvednumberofnodes+2];
    
    double x0 = nodecoordinates->at(3*node0+0); double y0 = nodecoordinates->at(3*node0+1); double z0 = nodecoordinates->at(3*node0+2);
    double x1 = nodecoordinates->at(3*node1+0); double y1 = nodecoordinates->at(3*node1+1); double z1 = nodecoordinates->at(3*node1+2);
    double x2 = nodecoordinates->at(3*node2+0); double y2 = nodecoordinates->at(3*node2+1); double z2 = nodecoordinates->at(3*node2+2);
    
    std::vector<double> a = {x1-x0, y1-y0, z1-z0};
    std::vector<double> b = {x2-x1, y2-y1, z2-z1};
    
    std::vector<double> crossprod = {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
    
    return crossprod;
}

int elements::getdimension(void)
{
    int maxi = 0;
    for (int i = 0; i < 8; i++)
    {
        if (count(i) > 0)
            maxi = i;
    }
    element elem(maxi);
    return elem.getelementdimension();
}

void elements::printnumber(void)
{
    std::cout << std::endl;
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        element myelement(elementtypenumber);
        std::cout << "Number of " << std::left << std::setw(12) << myelement.gettypenameconjugation(2) << " " << count(elementtypenumber) << std::endl;
    }
    std::cout << std::endl;
}

void elements::printsubelements(void)
{
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        element myelement(elementtypenumber, mycurvatureorder);

        if (subelementsinelements[elementtypenumber][0].size() != 0)
        {
            std::cout << "Points in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][0].size()/myelement.countcurvednodes(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.countcurvednodes(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][0][i*myelement.countcurvednodes()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (subelementsinelements[elementtypenumber][1].size() != 0)
        {
            std::cout << "Lines in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][1].size()/myelement.countedges(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.countedges(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][1][i*myelement.countedges()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (subelementsinelements[elementtypenumber][2].size() != 0)
        {
            std::cout << "Triangles in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][2].size()/myelement.counttriangularfaces(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.counttriangularfaces(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][2][i*myelement.counttriangularfaces()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        if (subelementsinelements[elementtypenumber][3].size() != 0)
        {
            std::cout << "Quadrangles in " << myelement.gettypenameconjugation(2) << ":" << std::endl;
            for (int i = 0; i < subelementsinelements[elementtypenumber][3].size()/myelement.countquadrangularfaces(); i++)
            {
                std::cout << std::left << std::setw(8) << i << " ";
                for (int j = 0; j < myelement.countquadrangularfaces(); j++)
                    std::cout << std::left << std::setw(6) << subelementsinelements[elementtypenumber][3][i*myelement.countquadrangularfaces()+j] << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}

void elements::printtotalorientations(void)
{
    std::cout << std::endl;
    for (int elementtypenumber = 0; elementtypenumber <= 7; elementtypenumber++)
    {
        if (totalorientations[elementtypenumber].size() == 0)
            continue;
        
        element myelement(elementtypenumber);
        std::cout << "Total orientations for all " << myelement.gettypenameconjugation(2) << ":" << std::endl;
        
        for (int i = 0; i < totalorientations[elementtypenumber].size(); i++)
            std::cout << std::left << std::setw(12) << i << std::left << std::setw(12) << totalorientations[elementtypenumber][i] << std::endl;
    }
    std::cout << std::endl;
}

void elements::write(std::string filename, int elementtypenumber, std::vector<int> elementnumbers)
{
    int numelems = elementnumbers.size();
    
    element myelem(elementtypenumber, mycurvatureorder);
    int ncn = myelem.countcurvednodes();
    
    densematrix xcoords(numelems, ncn);
    densematrix ycoords(numelems, ncn);
    densematrix zcoords(numelems, ncn);
    
    densematrix vals(numelems, ncn);

    double* xptr = xcoords.getvalues();
    double* yptr = ycoords.getvalues();
    double* zptr = zcoords.getvalues();
    double* vptr = vals.getvalues();
    
    for (int e = 0; e < numelems; e++)
    {
        std::vector<double> xc = getnodecoordinates(elementtypenumber, elementnumbers[e], 0);
        std::vector<double> yc = getnodecoordinates(elementtypenumber, elementnumbers[e], 1);
        std::vector<double> zc = getnodecoordinates(elementtypenumber, elementnumbers[e], 2);
        
        for (int n = 0; n < ncn; n++)
        {
            xptr[e*ncn+n] = xc[n];
            yptr[e*ncn+n] = yc[n];
            zptr[e*ncn+n] = zc[n];
            vptr[e*ncn+n] = elementnumbers[e];
        }
    }
    
    iodata datatowrite(mycurvatureorder, mycurvatureorder, true, {});
    
    datatowrite.addcoordinates(elementtypenumber, xcoords, ycoords, zcoords);
    datatowrite.adddata(elementtypenumber, {vals});
    
    iointerface::writetofile(filename, datatowrite);
}

std::vector<double> elements::computebarycenters(int elementtypenumber)
{
    std::vector<double>* nodecoordinates = mynodes->getcoordinates();

    // The barycenter of point elements is the nodes coordinates:
    if (elementtypenumber == 0)
        return *nodecoordinates;
        
    element myelement(elementtypenumber, mycurvatureorder);
    int nn = myelement.countnodes();
    int ncn = myelement.countcurvednodes();
    
    double invnn = 1.0/nn;
    
    // Preallocate the barycenter coordinates vector:
    std::vector<double> barycentercoordinates(3 * count(elementtypenumber),0);
    // Compute the barycenters on the straight element:
    int numel = count(elementtypenumber);
    for (int elem = 0; elem < numel; elem++)
    {
        for (int node = 0; node < nn; node++)
        {
            barycentercoordinates[3*elem+0] += nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elem*ncn+node]+0);
            barycentercoordinates[3*elem+1] += nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elem*ncn+node]+1);
            barycentercoordinates[3*elem+2] += nodecoordinates->at(3*subelementsinelements[elementtypenumber][0][elem*ncn+node]+2);
        }
        barycentercoordinates[3*elem+0] *= invnn;
        barycentercoordinates[3*elem+1] *= invnn;
        barycentercoordinates[3*elem+2] *= invnn;
    }
    
    return barycentercoordinates;
}

std::vector<int> elements::removeduplicates(int elementtypenumber)
{
    // For point elements (i.e. to remove node duplicates):
    if (elementtypenumber == 0)
    {
        std::vector<int> renumberingvector = mynodes->removeduplicates();
        renumber(0, renumberingvector);
        return renumberingvector;
    }
    
    int numberofcurvednodes = numberofsubelementsineveryelement[elementtypenumber][0];
    int numberoflines = numberofsubelementsineveryelement[elementtypenumber][1];
    int numberoftriangles = numberofsubelementsineveryelement[elementtypenumber][2];
    int numberofquadrangles = numberofsubelementsineveryelement[elementtypenumber][3];
    
    std::vector<double> barycentercoordinates = computebarycenters(elementtypenumber);
    
    // 'elementrenumbering' will give the renumbering corresponding to removed duplicates:
    std::vector<int> elementrenumbering;
    int numberofnonduplicates = myalgorithm::removeduplicates(barycentercoordinates, elementrenumbering);
    
    for (int i = 0; i < elementrenumbering.size(); i++)
    {
        if (elementrenumbering[i] != i)
        {
            for (int j = 0; j < numberofcurvednodes; j++)
                subelementsinelements[elementtypenumber][0][elementrenumbering[i]*numberofcurvednodes+j] = subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j];
            // Watch out: 'subelementsinelements[1][1]' does not contain anything!
            if (subelementsinelements[elementtypenumber][1].size() != 0)
            {
                for (int j = 0; j < numberoflines; j++)
                    subelementsinelements[elementtypenumber][1][elementrenumbering[i]*numberoflines+j] = subelementsinelements[elementtypenumber][1][i*numberoflines+j];
            }
            if (subelementsinelements[elementtypenumber][2].size() != 0)
            {
                for (int j = 0; j < numberoftriangles; j++)
                    subelementsinelements[elementtypenumber][2][elementrenumbering[i]*numberoftriangles+j] = subelementsinelements[elementtypenumber][2][i*numberoftriangles+j];
            }
            if (subelementsinelements[elementtypenumber][3].size() != 0)
            {
                for (int j = 0; j < numberofquadrangles; j++)
                    subelementsinelements[elementtypenumber][3][elementrenumbering[i]*numberofquadrangles+j] = subelementsinelements[elementtypenumber][3][i*numberofquadrangles+j];
            }
        }
    }
    // Shrink to fit:
    subelementsinelements[elementtypenumber][0].resize(numberofnonduplicates*numberofcurvednodes);
    if (subelementsinelements[elementtypenumber][1].size() != 0)
        subelementsinelements[elementtypenumber][1].resize(numberofnonduplicates*numberoflines);
    if (subelementsinelements[elementtypenumber][2].size() != 0)
        subelementsinelements[elementtypenumber][2].resize(numberofnonduplicates*numberoftriangles);
    if (subelementsinelements[elementtypenumber][3].size() != 0)
        subelementsinelements[elementtypenumber][3].resize(numberofnonduplicates*numberofquadrangles);

    renumber(elementtypenumber, elementrenumbering);
    
    return elementrenumbering;
}

void elements::renumber(int elementtypenumber, std::vector<int>& renumberingvector)
{   
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        switch (elementtypenumber)
        {
            // Renumber the point elements (i.e. the nodes):
            case 0:
                for (int i = 0; i < subelementsinelements[typenum][0].size(); i++)
                    subelementsinelements[typenum][0][i] = renumberingvector[subelementsinelements[typenum][0][i]];
                break;
            // Renumber the line elements:
            case 1:
                for (int i = 0; i < subelementsinelements[typenum][1].size(); i++)
                    subelementsinelements[typenum][1][i] = renumberingvector[subelementsinelements[typenum][1][i]];
                break;
            // Renumber the triangle elements:
            case 2:
                for (int i = 0; i < subelementsinelements[typenum][2].size(); i++)
                    subelementsinelements[typenum][2][i] = renumberingvector[subelementsinelements[typenum][2][i]];
                break;
            // Renumber the quadrangle elements:
            case 3:
                for (int i = 0; i < subelementsinelements[typenum][3].size(); i++)
                    subelementsinelements[typenum][3][i] = renumberingvector[subelementsinelements[typenum][3][i]];
                break;
        }
    }
}

void elements::reorder(int elementtypenumber, std::vector<int> &elementreordering)
{
    // For point elements (i.e. nodes):
    if (elementtypenumber == 0)
        mynodes->reorder(elementreordering);
    else
    {
        int numberofnodes = numberofsubelementsineveryelement[elementtypenumber][0];
        int numberoflines = numberofsubelementsineveryelement[elementtypenumber][1];
        int numberoftriangles = numberofsubelementsineveryelement[elementtypenumber][2];
        int numberofquadrangles = numberofsubelementsineveryelement[elementtypenumber][3];
        
        std::vector<int> pointsinelementspart = subelementsinelements[elementtypenumber][0];
        for (int i = 0; i < subelementsinelements[elementtypenumber][0].size()/numberofnodes; i++)
        {
            for (int j = 0; j < numberofnodes; j++)
                subelementsinelements[elementtypenumber][0][numberofnodes*i+j] = pointsinelementspart[numberofnodes*elementreordering[i]+j]; 
        }
        std::vector<int> linesinelementspart = subelementsinelements[elementtypenumber][1];
        for (int i = 0; i < subelementsinelements[elementtypenumber][1].size()/numberoflines; i++)
        {
            for (int j = 0; j < numberoflines; j++)
                subelementsinelements[elementtypenumber][1][numberoflines*i+j] = linesinelementspart[numberoflines*elementreordering[i]+j]; 
        }
        if (numberoftriangles != 0)
        {
            std::vector<int> trianglesinelementspart = subelementsinelements[elementtypenumber][2];
            for (int i = 0; i < subelementsinelements[elementtypenumber][2].size()/numberoftriangles; i++)
            {
                for (int j = 0; j < numberoftriangles; j++)
                    subelementsinelements[elementtypenumber][2][numberoftriangles*i+j] = trianglesinelementspart[numberoftriangles*elementreordering[i]+j]; 
            }
        }
        if (numberofquadrangles != 0)
        {
            std::vector<int> quadranglesinelementspart = subelementsinelements[elementtypenumber][3];
            for (int i = 0; i < subelementsinelements[elementtypenumber][3].size()/numberofquadrangles; i++)
            {
                for (int j = 0; j < numberofquadrangles; j++)
                    subelementsinelements[elementtypenumber][3][numberofquadrangles*i+j] = quadranglesinelementspart[numberofquadrangles*elementreordering[i]+j]; 
            }
        }
    }
    
    
    std::vector<int> indisjointregionpart = indisjointregion[elementtypenumber];
    for (int i = 0; i < indisjointregion[elementtypenumber].size(); i++)
        indisjointregion[elementtypenumber][i] = indisjointregionpart[elementreordering[i]];

    std::vector<int> totalorientationspart = totalorientations[elementtypenumber];
    for (int i = 0; i < totalorientations[elementtypenumber].size(); i++)
        totalorientations[elementtypenumber][i] = totalorientationspart[elementreordering[i]];

    std::vector<double> barycenterspart = barycenters[elementtypenumber];
    for (int i = 0; i < barycenters[elementtypenumber].size()/3; i++)
    {
        barycenters[elementtypenumber][3*i+0] = barycenterspart[3*elementreordering[i]+0];
        barycenters[elementtypenumber][3*i+1] = barycenterspart[3*elementreordering[i]+1];
        barycenters[elementtypenumber][3*i+2] = barycenterspart[3*elementreordering[i]+2];
    }
    
    std::vector<double> sphereradiuspart = sphereradius[elementtypenumber];
    for (int i = 0; i < sphereradius[elementtypenumber].size(); i++)
        sphereradius[elementtypenumber][i] = sphereradiuspart[elementreordering[i]];
    
    std::vector<double> boxdimensionspart = boxdimensions[elementtypenumber];
    for (int i = 0; i < boxdimensions[elementtypenumber].size()/3; i++)
    {
        boxdimensions[elementtypenumber][3*i+0] = boxdimensionspart[3*elementreordering[i]+0];
        boxdimensions[elementtypenumber][3*i+1] = boxdimensionspart[3*elementreordering[i]+1];
        boxdimensions[elementtypenumber][3*i+2] = boxdimensionspart[3*elementreordering[i]+2];
    }
    
    adressedgesatnodes = {};
    edgesatnodes = {};
    
    adresscellsatedges = {};
    cellsatedges = {};
}

void elements::explode(void)
{                            
    // Add all new elements. Loop on all elements with increasing dimension 
    // to avoid defining too many duplicates.
    // Skip lines (type 1) since there is nothing to add for lines.
    for (int elementtypenumber = 2; elementtypenumber <= 7; elementtypenumber++)
    {
        element myelement(elementtypenumber, mycurvatureorder);
        int curvednumberofnodes = myelement.countcurvednodes();
        int numberofedges = myelement.countedges();
        int numberoftriangularfaces = myelement.counttriangularfaces();
        int numberofquadrangularfaces = myelement.countquadrangularfaces();
        int elementdimension = myelement.getelementdimension();
        
        std::vector<int> facesdefinitionsbasedonedges = myelement.getfacesdefinitionsbasedonedges();

        // Skip if there are no elements of the current type:
        if (count(elementtypenumber) == 0)
            continue;

        std::vector<int> currentnodes(curvednumberofnodes,0);
        
        // Get all line/triangle/quadrangle definitions:
        std::vector<int> consecutives = myalgorithm::getequallyspaced(0,1,curvednumberofnodes);
        myelement.setnodes(consecutives);
        // Extract line node indexes:
        std::vector<std::vector<int>> indexesinlines(numberofedges);
        for (int i = 0; i < numberofedges; i++)
            indexesinlines[i] = myelement.getnodesinline(i);
        int numnodesinline = indexesinlines[0].size();
        std::vector<int> nodesinline(numnodesinline);
        // Extract triangular face node indexes:
        std::vector<std::vector<int>> indexesintris(numberoftriangularfaces);
        for (int i = 0; i < numberoftriangularfaces; i++)
            indexesintris[i] = myelement.getnodesintriangle(i);
        int numnodesintri = 0;
        if (numberoftriangularfaces > 0)
            numnodesintri = indexesintris[0].size();
        std::vector<int> nodesintriangle(numnodesintri);
        // Extract quadrangular face node indexes:
        std::vector<std::vector<int>> indexesinquads(numberofquadrangularfaces);
        for (int i = 0; i < numberofquadrangularfaces; i++)
            indexesinquads[i] = myelement.getnodesinquadrangle(i);
        int numnodesinquad = 0;
        if (numberofquadrangularfaces > 0)
            numnodesinquad = indexesinquads[0].size();
        std::vector<int> nodesinquadrangle(numnodesinquad);

        // Loop on all elements of the current type:
        for (int elem = 0; elem < count(elementtypenumber); elem++)
        {
            // Define the nodes of the element object:
            for (int i = 0; i < curvednumberofnodes; i++)
                currentnodes[i] = subelementsinelements[elementtypenumber][0][elem*curvednumberofnodes+i];

            // Add all line subelements - only for 2D and 3D elements:
            int firstnewlinenumber = count(1);
            for (int line = 0; line < numberofedges; line++)
            {
                for (int i = 0; i < numnodesinline; i++)
                    nodesinline[i] = currentnodes[indexesinlines[line][i]];
                    
                int newlinenumber = add(1, mycurvatureorder, nodesinline);
                subelementsinelements[elementtypenumber][1].push_back(newlinenumber);
            }
            // Add all subsurfaces - only for 3D elements:
            if (elementdimension == 3)
            {
                for (int triangle = 0; triangle < numberoftriangularfaces; triangle++)
                {
                    for (int i = 0; i < numnodesintri; i++)
                        nodesintriangle[i] = currentnodes[indexesintris[triangle][i]];
                    
                    int newtrianglenumber = add(2, mycurvatureorder, nodesintriangle);
                    subelementsinelements[elementtypenumber][2].push_back(newtrianglenumber);
                    // Also link the triangle to its lines. All lines have already been defined before above.
                    for (int triangleline = 0; triangleline < 3; triangleline++)
                        subelementsinelements[2][1].push_back(firstnewlinenumber + std::abs(facesdefinitionsbasedonedges[3*triangle+triangleline])-1);
                }    
                for (int quadrangle = 0; quadrangle < numberofquadrangularfaces; quadrangle++)
                {
                    for (int i = 0; i < numnodesinquad; i++)
                        nodesinquadrangle[i] = currentnodes[indexesinquads[quadrangle][i]];
                        
                    int newquadranglenumber = add(3, mycurvatureorder, nodesinquadrangle);
                    subelementsinelements[elementtypenumber][3].push_back(newquadranglenumber);
                    // Also link the quadrangle to its lines (defined before):
                    for (int quadrangleline = 0; quadrangleline < 4; quadrangleline++)
                        subelementsinelements[3][1].push_back(firstnewlinenumber + std::abs(facesdefinitionsbasedonedges[3*numberoftriangularfaces+4*quadrangle+quadrangleline])-1);
                }    
            }
        }
    }
}

void elements::definedisjointregions(void)
{
    
    int numberofphysicalregions = myphysicalregions->count();
    
    // Define 'isinphysicalregion[elementtypenumber]' to hold a vector
    // whose (elem*numphysreg+i)th entry is true if the element 'elem' 
    // is included in the ith physical region.
    std::vector<std::vector<bool>> isinphysicalregion(8);
    for (int typenum = 0; typenum <= 7; typenum++)
        isinphysicalregion[typenum] = std::vector<bool>(count(typenum)*numberofphysicalregions);
    
    // Write the physical regions to the elements. 
    for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
    {
        physicalregion* currentphysicalregion = myphysicalregions->getatindex(physregindex);
        std::vector<std::vector<int>>* elementsinphysicalregion = currentphysicalregion->getelementlist();
        
        for (int typenum = 0; typenum <= 7; typenum++)
        {
            // Iterate on all elements of the given type:
            for (int i = 0; i < (*elementsinphysicalregion)[typenum].size(); i++)
                isinphysicalregion[typenum][(*elementsinphysicalregion)[typenum][i] * numberofphysicalregions + physregindex] = true;
        }
    }
        
    
    // Propagate the physical region memberships from elements to the subelements.
    for (int typenum = 1; typenum <= 7; typenum++)
    {
        int numberofcurvednodes = numberofsubelementsineveryelement[typenum][0];
        int numberoflines = numberofsubelementsineveryelement[typenum][1];
        int numberoftriangles = numberofsubelementsineveryelement[typenum][2];
        int numberofquadrangles = numberofsubelementsineveryelement[typenum][3];

        for (int elem = 0; elem < subelementsinelements[typenum][0].size()/numberofcurvednodes; elem++)
        {
            for (int i = 0; i < numberofcurvednodes; i++)
            {
                int subelem = subelementsinelements[typenum][0][elem*numberofcurvednodes+i];
                for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                {
                    if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                        isinphysicalregion[0][subelem*numberofphysicalregions+physregindex] = true;
                }
            }
        }
        
        for (int elem = 0; elem < subelementsinelements[typenum][1].size()/numberoflines; elem++)
        {
            for (int i = 0; i < numberoflines; i++)
            {
                int subelem = subelementsinelements[typenum][1][elem*numberoflines+i];
                for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                {
                    if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                        isinphysicalregion[1][subelem*numberofphysicalregions+physregindex] = true;
                }
            }
        }

        if (numberoftriangles != 0)
        {
            for (int elem = 0; elem < subelementsinelements[typenum][2].size()/numberoftriangles; elem++)
            {
                for (int i = 0; i < numberoftriangles; i++)
                {
                    int subelem = subelementsinelements[typenum][2][elem*numberoftriangles+i];
                    for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                    {
                        if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                            isinphysicalregion[2][subelem*numberofphysicalregions+physregindex] = true;
                    }
                }
            }
        }
        
        if (numberofquadrangles != 0)
        {
            for (int elem = 0; elem < subelementsinelements[typenum][3].size()/numberofquadrangles; elem++)
            {
                for (int i = 0; i < numberofquadrangles; i++)
                {
                    int subelem = subelementsinelements[typenum][3][elem*numberofquadrangles+i];
                    for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                    {
                        if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                            isinphysicalregion[3][subelem*numberofphysicalregions+physregindex] = true;
                    }
                }
            }
        }
    }
        

    // Now define the disjoint regions based on 'isinphysicalregion':
    indisjointregion.resize(8);
    for (int typenum = 0; typenum <= 7; typenum++)
        indisjointregion[typenum].resize(count(typenum));
    // We need to know if a node is a corner node or not:
    std::vector<bool> isnodeacornernode = iscornernode();
    // Define the disjoint regions and fill in 'indisjointregion':
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        if (count(typenum) == 0)
            continue;

        std::vector<bool> temp(numberofphysicalregions);    

        for (int elem = 0; elem < count(typenum); elem++)
        {   
            // Only corner nodes matter for the disjoint regions 
            // (no dof will ever be associated to curvature nodes).
            if (typenum != 0 || isnodeacornernode[elem])
            {
                for (int i = 0; i < numberofphysicalregions; i++)
                    temp[i] = isinphysicalregion[typenum][elem*numberofphysicalregions+i];
                indisjointregion[typenum][elem] = mydisjointregions->add(typenum, temp);
            }
            else
                indisjointregion[typenum][elem] = -1;
        }
    }
}

void elements::reorderbydisjointregions(void)
{
    std::vector<std::vector<int>> elementrenumbering;
    reorderbydisjointregions(elementrenumbering);
}

void elements::reorderbydisjointregions(std::vector<std::vector<int>>& elementrenumbering)
{
    elementrenumbering.resize(8);
    
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        std::vector<int> elementreordering;
        myalgorithm::stablesort(indisjointregion[typenum], elementreordering);
        
        elementrenumbering[typenum] = std::vector<int>(count(typenum));
        for (int i = 0; i < count(typenum); i++)
            elementrenumbering[typenum][elementreordering[i]] = i;
        
        reorder(typenum, elementreordering);
        renumber(typenum, elementrenumbering[typenum]);

        for (int physregindex = 0; physregindex < myphysicalregions->count(); physregindex++)
        {
            physicalregion* currentphysicalregion = myphysicalregions->getatindex(physregindex);
            currentphysicalregion->renumberelements(typenum, elementrenumbering[typenum]);
        }
    }
}

void elements::definedisjointregionsranges(void)
{
    for (int typenum = 0; typenum <= 7; typenum++)
    {
        for (int i = 0; i < indisjointregion[typenum].size(); i++)
        {
            if (i == 0)
                mydisjointregions->setrangebegin(indisjointregion[typenum][i], 0);
            if (i != 0 && indisjointregion[typenum][i] != indisjointregion[typenum][i-1])
            {
                mydisjointregions->setrangeend(indisjointregion[typenum][i-1], i-1);
                mydisjointregions->setrangebegin(indisjointregion[typenum][i], i);
            }
            if (i == indisjointregion[typenum].size()-1)
                mydisjointregions->setrangeend(indisjointregion[typenum][i], indisjointregion[typenum].size()-1);
        }
    }
}

std::vector<bool> elements::iscornernode(void)
{
    std::vector<bool> output(count(0), false);
    for (int typenum = 1; typenum <= 7; typenum++)
    {
        element curvedelement(typenum, getcurvatureorder());
        int numcurvednodes = curvedelement.countcurvednodes();
        int numcornernodes = curvedelement.countnodes();

        for (int i = 0; i < count(typenum); i++)
        {    
            // The first 'numcornernodes' are corner nodes, the remaining ones not.
            for (int node = 0; node < numcornernodes; node++)
                output[subelementsinelements[typenum][0][i*numcurvednodes+node]] = true;
        }
    }
    return output;
}

void elements::orient(void)
{
    // Loop on all element types except the point element (type 0):
    for (int elementtypenumber = 1; elementtypenumber <= 7; elementtypenumber++)
    {    
        if (subelementsinelements[elementtypenumber][0].size() == 0)
            continue;
        
        element myelement(elementtypenumber, mycurvatureorder);
        
        int numberofnodes = myelement.countnodes();
        int numberofcurvednodes = myelement.countcurvednodes();
        int numelemofcurrenttype = count(elementtypenumber);
        
        // Loop on all elements:
        std::vector<int> cornernodes(numelemofcurrenttype*numberofnodes);
        for (int i = 0; i < numelemofcurrenttype; i++)
        {
            for (int j = 0; j < numberofnodes; j++)
                cornernodes[i*numberofnodes+j] = subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j];
        }
        totalorientations[elementtypenumber] = orientation::gettotalorientation(elementtypenumber, cornernodes); 
    }
}

elements elements::copy(nodes* nds, physicalregions* prs, disjointregions* drs)
{
    elements out;
    
    out = *this;
    
    out.mynodes = nds;
    out.myphysicalregions = prs;
    out.mydisjointregions = drs;
    
    return out;
}

void elements::toptracker(std::shared_ptr<ptracker> originpt, std::shared_ptr<ptracker> targetpt)
{
    std::vector<std::vector<int>> renumbering;
    std::vector<std::vector<int>> reordering(8, std::vector<int>(0));
    originpt->getrenumbering(targetpt, renumbering);
    
    for (int i = 0; i < 8; i++)
    {
        if (renumbering[i].size() == 0)
            continue;
    
        reordering[i] = myalgorithm::getreordering(renumbering[i]);    
        // Renumber and reorder the elements in all containers:
        renumber(i, renumbering[i]);
        reorder(i, reordering[i]);  
    }
    
    targetpt->getindisjointregions(indisjointregion);
}

