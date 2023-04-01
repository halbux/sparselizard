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
        logs log;
        log.msg() << "Error in 'elements' object: the mesh can contain only a single curvature order." << std::endl;
        log.error();
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
        logs log;
        log.msg() << "Error in 'elements' object: returned a negative disjoint region number" << std::endl;
        log.msg() << "Possible reason: too bad quality elements obtained when h-adapting a curved mesh (some identical edges might not have been merged and therefore AMR has failed)" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int elements::gettotalorientation(int elementtypenumber, int elementnumber)
{
    if (elementtypenumber != 0)
        return totalorientations[elementtypenumber][elementnumber];
    else
        return 0;
}

std::vector<bool> elements::isflipped(int subelementtypenumber, std::vector<int>& subelementnumbers, int elementtypenumber, std::vector<int>& elementnumbers)
{
    int numelems = elementnumbers.size();
    std::vector<bool> output(numelems, false);

    // Points cannot be flipped:
    if (subelementtypenumber == 0)
        return output;

    element parent(elementtypenumber, mycurvatureorder);
    element sub(subelementtypenumber, mycurvatureorder);
    
    int numsubsinparent = parent.counttype(subelementtypenumber);
    int numcurvednodesinparent = parent.countcurvednodes();
    int numtrisinparent = parent.counttriangularfaces();
    
    int numnodesinsub = sub.countnodes();
    int numcurvednodesinsub = sub.countcurvednodes();

    std::vector<int> cornernodesinparent;
    if (subelementtypenumber == 1)
        cornernodesinparent = parent.getedgesdefinitionsbasedonnodes();
    else
        cornernodesinparent = parent.getfacesdefinitionsbasedonnodes();
    int trishift = 0;
    if (subelementtypenumber == 3)
        trishift = 3*numtrisinparent;
    
    for (int e = 0; e < numelems; e++)
    {
        int curparent = elementnumbers[e];
        int cursub = subelementnumbers[e];
        
        // Find the subelement index in the parent:
        int subindex = -1;
        for (int i = 0; i < numsubsinparent; i++)
        {
            int subnum = subelementsinelements[elementtypenumber][subelementtypenumber][curparent*numsubsinparent+i]; // subs are never points thanks to above check
            if (cursub == subnum)
                subindex = i;
        }
        if (subindex == -1)
        {
            logs log;
            log.msg() << "Error in 'elements' object: in 'isflipped' could not find " << sub.gettypename() <<  " " << cursub << " in " << parent.gettypename() << " " << curparent << std::endl;
            log.error();
        }
        
        // Extract the corner nodes of the sub:
        std::vector<int> subcorners(numnodesinsub);
        for (int i = 0; i < numnodesinsub; i++)
            subcorners[i] = subelementsinelements[subelementtypenumber][0][cursub*numcurvednodesinsub+i];
        
        // Extract the corner nodes of the subelement in the parent:
        int firstcorner = trishift+numnodesinsub*subindex;
        
        std::vector<int> subcornersinparent(numnodesinsub);
        for (int i = 0; i < numnodesinsub; i++)
            subcornersinparent[i] = subelementsinelements[elementtypenumber][0][curparent*numcurvednodesinparent + cornernodesinparent[firstcorner+i]];

        // Compare orientation between sub and sub in parent:
        output[e] = gentools::isflipped(subcorners, subcornersinparent);
    }

    return output;
}

int elements::istypeinelementlists(int elementtypenumber, std::vector<std::vector<std::vector<int>>*> elementlists, std::vector<bool>& isinelementlists, bool considercurvaturenodes)
{
    isinelementlists = std::vector<bool>(count(elementtypenumber), false);

    int numinlists = 0;
    for (int i = 0; i < 8; i++)
    {
        element el(i, mycurvatureorder);
        int ns;
        if (elementtypenumber == 0 && considercurvaturenodes)
            ns = el.countcurvednodes();
        else
            ns = el.counttype(elementtypenumber);
            
        if (ns == 0)
            continue;

        for (int l = 0; l < elementlists.size(); l++)
        {
            if (elementlists[l] != NULL)
            {
                for (int j = 0; j < elementlists[l]->at(i).size(); j++)
                {
                    int curelem = elementlists[l]->at(i)[j];
                    for (int k = 0; k < ns; k++)
                    {
                        int cursub = getsubelement(elementtypenumber,i,curelem,k);
                        if (isinelementlists[cursub] == false)
                        {
                            isinelementlists[cursub] = true;
                            numinlists++;
                        }
                    }
                }
            }
        }
    }

    return numinlists;
}

int elements::istypeindisjointregions(int elementtypenumber, std::vector<bool> isdisjregselected, std::vector<bool>& isinelementlists, bool considercurvaturenodes)
{
    isinelementlists = std::vector<bool>(count(elementtypenumber), false);

    int numinlists = 0;
    for (int d = 0; d < isdisjregselected.size(); d++)
    {
        if (isdisjregselected[d] == false)
            continue;
            
        int tn = mydisjointregions->getelementtypenumber(d);
        int rb = mydisjointregions->getrangebegin(d);
        int ne = mydisjointregions->countelements(d);
    
        element el(tn, mycurvatureorder);
        int ns;
        if (elementtypenumber == 0 && considercurvaturenodes)
            ns = el.countcurvednodes();
        else
            ns = el.counttype(elementtypenumber);
            
        if (ns == 0)
            continue;

        for (int e = 0; e < ne; e++)
        {
            for (int k = 0; k < ns; k++)
            {
                int cursub = getsubelement(elementtypenumber,tn,rb+e,k);
                if (isinelementlists[cursub] == false)
                {
                    isinelementlists[cursub] = true;
                    numinlists++;
                }
            }
        }
    }

    return numinlists;
}

int elements::count(int elementtypenumber)
{
    if (elementtypenumber == 0)
        return mynodes->count();
    
    if (subelementsinelements[elementtypenumber][0].size() == 0)
        return 0;
    
    return subelementsinelements[elementtypenumber][0].size()/numberofsubelementsineveryelement[elementtypenumber][0];
}

int elements::countcells(int elementtypenumber)
{
    element el(elementtypenumber);
    if (el.getelementdimension() != getdimension())
        return 0;
    else
        return count(elementtypenumber);
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

void elements::populatecellsattype(int subtype, std::vector<int>& act, std::vector<int>& ct)
{
    int celldim = getdimension();

    // Get the number of subs:
    int numsubs = count(subtype);
    // Get the number of cells of each type and the number of subs in each cell:
    std::vector<int> numcells(8,0);
    std::vector<int> ns(8);
    int prealloc = 0;
    for (int i = 0; i < 8; i++)
    {
        element myelem(i);
        if (myelem.getelementdimension() == celldim)
        {
            numcells[i] = count(i);
            ns[i] = myelem.counttype(subtype);
            prealloc += numcells[i]*ns[i];   
        }
    }


    // Preallocate vectors:
    act = std::vector<int>(numsubs,0);
    ct = std::vector<int>(2*prealloc);


    // Vector to count the number of cells touching every sub:
    std::vector<int> numcellsonsub(numsubs,0);

    // Loop on all cells:
    for (int i = 0; i < 8; i++)
    {
        for (int c = 0; c < numcells[i]; c++)
        {
            for (int s = 0; s < ns[i]; s++)
                numcellsonsub[getsubelement(subtype,i,c,s)]++;
        }
    }

    for (int s = 1; s < numsubs; s++)
        act[s] = act[s-1] + 2*numcellsonsub[s-1];

    std::vector<int> currentcellinsub(numsubs,0); 
    for (int i = 0; i < 8; i++)
    {
        for (int c = 0; c < numcells[i]; c++)
        {
            for (int s = 0; s < ns[i]; s++)
            {
                int cursub = getsubelement(subtype,i,c,s);
                ct[act[cursub]+2*currentcellinsub[cursub]+0] = i;
                ct[act[cursub]+2*currentcellinsub[cursub]+1] = c;
                
                currentcellinsub[cursub]++;
            }
        }
    }
}

int elements::countcellsontype(int subtype, int subnumber)
{
    if (adresscellsattype[subtype].size() == 0)
        populatecellsattype(subtype, adresscellsattype[subtype], cellsattype[subtype]);

    if (subnumber+1 < adresscellsattype[subtype].size())
        return (adresscellsattype[subtype][subnumber+1]-adresscellsattype[subtype][subnumber])/2;
    else
        return (cellsattype[subtype].size()-adresscellsattype[subtype][subnumber])/2;
}

std::vector<int> elements::getcellsontype(int subtype, int subnumber)
{
    int numcellsatsub = countcellsontype(subtype, subnumber);

    std::vector<int> output(2*numcellsatsub);
    for (int i = 0; i < 2*numcellsatsub; i++)
        output[i] = cellsattype[subtype][adresscellsattype[subtype][subnumber]+i];

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

void elements::getbarycenters(std::vector<std::vector<int>>* elementlist, std::vector<double>& barys)
{
    int numelems = 0;
    for (int i = 0; i < 8; i++)
        numelems += elementlist->at(i).size();

    barys = std::vector<double>(3*numelems);

    int index = 0;
    for (int i = 0; i < 8; i++)
    {
        if (elementlist->at(i).size() == 0)
            continue;

        double* bcs = getbarycenters(i)->data();
        for (int j = 0; j < elementlist->at(i).size(); j++)
        {
            int elem = elementlist->at(i)[j];

            barys[3*index+0] = bcs[3*elem+0];
            barys[3*index+1] = bcs[3*elem+1];
            barys[3*index+2] = bcs[3*elem+2];

            index++;
        }
    }
}

void elements::getbarycenters(int elementtypenumber, std::vector<int>& elementlist, std::vector<double>& barys)
{
    barys = std::vector<double>(3*elementlist.size());
    getbarycenters(elementtypenumber, elementlist, barys.data());
}

void elements::getbarycenters(int elementtypenumber, std::vector<int>& elementlist, double* barys)
{
    int numelems = elementlist.size();
    
    if (numelems == 0)
        return;
    
    double* bcs = getbarycenters(elementtypenumber)->data();

    for (int i = 0; i < numelems; i++)
    {
        int elem = elementlist[i];

        barys[3*i+0] = bcs[3*elem+0];
        barys[3*i+1] = bcs[3*elem+1];
        barys[3*i+2] = bcs[3*elem+2];
    }
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

int elements::countinterfaceelementtypes(void)
{
    std::vector<int> countfordim = {-1,1,2,4};
    return countfordim[getdimension()];
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

void elements::write(std::string filename, int elementtypenumber, std::vector<int> elementnumbers, std::vector<int> elementvalues)
{
    int numelems = elementnumbers.size();
    if (numelems == 0)
        return;
    
    element myelem(elementtypenumber, mycurvatureorder);
    int ncn = myelem.countcurvednodes();
    
    densemat xcoords(numelems, ncn);
    densemat ycoords(numelems, ncn);
    densemat zcoords(numelems, ncn);
    
    densemat vals(numelems, ncn);

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
            vptr[e*ncn+n] = elementvalues[e];
        }
    }
    
    iodata datatowrite(mycurvatureorder, mycurvatureorder, true, {});
    
    datatowrite.addcoordinates(elementtypenumber, xcoords, ycoords, zcoords);
    datatowrite.adddata(elementtypenumber, {vals});
    
    iointerface::writetofile(filename, datatowrite);
}

void elements::writeedgedirection(int physreg, std::string filename)
{
    std::vector<int> ders = myphysicalregions->get(physreg)->getdisjointregions(1);

    int numedgesinpr = 0;
    for (int i = 0; i < ders.size(); i++)
        numedgesinpr += mydisjointregions->countelements(ders[i]);
    
    densemat xcoords(numedgesinpr, 1);
    densemat ycoords(numedgesinpr, 1);
    densemat zcoords(numedgesinpr, 1);
    
    densemat xvals(numedgesinpr, 1);
    densemat yvals(numedgesinpr, 1);
    densemat zvals(numedgesinpr, 1);
    
    double* xptr = xcoords.getvalues();
    double* yptr = ycoords.getvalues();
    double* zptr = zcoords.getvalues();
    
    double* vxptr = xvals.getvalues();
    double* vyptr = yvals.getvalues();
    double* vzptr = zvals.getvalues();
    
    double* nodecoords = mynodes->getcoordinates()->data();
    
    int index = 0;
    for (int i = 0; i < ders.size(); i++)
    {
        int rb = mydisjointregions->getrangebegin(ders[i]);
        int ne = mydisjointregions->countelements(ders[i]);
        
        for (int j = 0; j < ne; j++)
        {
            int firstnode = getsubelement(0, 1, rb+j, 0);
            int lastnode = getsubelement(0, 1, rb+j, 1);
        
            double xf = nodecoords[3*firstnode+0];
            double yf = nodecoords[3*firstnode+1];
            double zf = nodecoords[3*firstnode+2];
            
            double xl = nodecoords[3*lastnode+0];
            double yl = nodecoords[3*lastnode+1];
            double zl = nodecoords[3*lastnode+2];
        
            // Uncurved edge barycenter:
            xptr[index] = 0.5*(xf+xl);
            yptr[index] = 0.5*(yf+yl);
            zptr[index] = 0.5*(zf+zl);

            vxptr[index] = xl-xf;
            vyptr[index] = yl-yf;
            vzptr[index] = zl-zf;
            
            double nrm = std::sqrt(vxptr[index]*vxptr[index] + vyptr[index]*vyptr[index] + vzptr[index]*vzptr[index]);
            
            vxptr[index] /= nrm;
            vyptr[index] /= nrm;
            vzptr[index] /= nrm;
            
            index++;
        }
    }
    
    iodata datatowrite(mycurvatureorder, mycurvatureorder, false, {});
    
    datatowrite.addcoordinates(0, xcoords, ycoords, zcoords);
    datatowrite.adddata(0, {xvals, yvals, zvals});
    
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
    
    element elobj(elementtypenumber, mycurvatureorder);
    
    int numberofcurvednodes = elobj.countcurvednodes();
    int numberoflines = elobj.countedges();
    int numberoftriangles = elobj.counttriangularfaces();
    int numberofquadrangles = elobj.countquadrangularfaces();
    
    std::vector<int> elementrenumbering;
    int numberofnonduplicates = gentools::removeduplicates(subelementsinelements[elementtypenumber][0], elementrenumbering, numberofcurvednodes);
    
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
    
    adresscellsattype = std::vector<std::vector<int>>(4, std::vector<int>(0));
    cellsattype = std::vector<std::vector<int>>(4, std::vector<int>(0));
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
        std::vector<int> consecutives = gentools::getequallyspaced(0,1,curvednumberofnodes);
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

void elements::follow(std::vector<std::vector<int>>* elementlist, int subtype, std::vector<int>& sublist, std::vector<std::vector<std::vector<int>>*> mustbeinelementlists)
{
    int highestdim = gentools::getmaxdim(elementlist);

    int numsubs = count(subtype);
    
    // To treat only once each subelement:
    std::vector<bool> issubdone(numsubs, false);
    // Preallocate to max possible size:
    sublist.resize(numsubs);
    
    std::vector<bool> issuballowed = {};
    if (mustbeinelementlists.size() > 0)
        istypeinelementlists(subtype, mustbeinelementlists, issuballowed, false);

    int index = 0;
    for (int i = 0; i < 8; i++)
    {
        int numelems = elementlist->at(i).size();
    
        element el(i);
        int ns = el.counttype(subtype); // no curvature nodes
        
        if (numelems == 0 || ns == 0 || el.getelementdimension() != highestdim)
            continue;

        for (int j = 0; j < numelems; j++)
        {
            int elem = elementlist->at(i)[j];
            for (int k = 0; k < ns; k++)
            {
                int cursub = getsubelement(subtype,i,elem,k);

                if (issubdone[cursub] == false && (mustbeinelementlists.size() == 0 || issuballowed[cursub]))
                {
                    sublist[index] = cursub;
                    issubdone[cursub] = true;

                    index++;
                }
            }
        }
    }

    sublist.resize(index);
}

void elements::definedisjointregions(void)
{
    long long int numberofphysicalregions = myphysicalregions->count();
    
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
            for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
            {
                if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                {
                    for (int i = 0; i < numberofcurvednodes; i++)
                    {
                        int subelem = subelementsinelements[typenum][0][elem*numberofcurvednodes+i];
                        isinphysicalregion[0][subelem*numberofphysicalregions+physregindex] = true;
                    }
                }
            }
        }
        
        for (int elem = 0; elem < subelementsinelements[typenum][1].size()/numberoflines; elem++)
        {
            for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
            {
                if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                {
                    for (int i = 0; i < numberoflines; i++)
                    {
                        int subelem = subelementsinelements[typenum][1][elem*numberoflines+i];
                        isinphysicalregion[1][subelem*numberofphysicalregions+physregindex] = true;
                    }
                }
            }
        }

        if (numberoftriangles != 0)
        {
            for (int elem = 0; elem < subelementsinelements[typenum][2].size()/numberoftriangles; elem++)
            {
                for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                {
                    if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                    {
                        for (int i = 0; i < numberoftriangles; i++)
                        {
                            int subelem = subelementsinelements[typenum][2][elem*numberoftriangles+i];
                            isinphysicalregion[2][subelem*numberofphysicalregions+physregindex] = true;
                        }
                    }
                }
            }
        }
        
        if (numberofquadrangles != 0)
        {
            for (int elem = 0; elem < subelementsinelements[typenum][3].size()/numberofquadrangles; elem++)
            {
                for (int physregindex = 0; physregindex < numberofphysicalregions; physregindex++)
                {
                    if (isinphysicalregion[typenum][elem*numberofphysicalregions+physregindex])
                    {
                        for (int i = 0; i < numberofquadrangles; i++)
                        {
                            int subelem = subelementsinelements[typenum][3][elem*numberofquadrangles+i];
                            isinphysicalregion[3][subelem*numberofphysicalregions+physregindex] = true;
                        }
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
        gentools::stablesort(indisjointregion[typenum], elementreordering);
        
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
    std::vector<bool> output(count(0), true);
    for (int typenum = 1; typenum <= 7; typenum++)
    {
        element curvedelement(typenum, getcurvatureorder());
        int numcurvednodes = curvedelement.countcurvednodes();
        int numcornernodes = curvedelement.countnodes();

        for (int i = 0; i < count(typenum); i++)
        {    
            // The first 'numcornernodes' are corner nodes, the remaining ones not.
            for (int node = numcornernodes; node < numcurvednodes; node++)
                output[subelementsinelements[typenum][0][i*numcurvednodes+node]] = false;
        }
    }
    return output;
}

void elements::orient(long long int* noderenumbering)
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
        std::vector<long long int> cornernodes(numelemofcurrenttype*numberofnodes);
        for (int i = 0; i < numelemofcurrenttype; i++)
        {
            for (int j = 0; j < numberofnodes; j++)
            {
                int curnode = subelementsinelements[elementtypenumber][0][i*numberofcurvednodes+j];
                if (noderenumbering == NULL)
                    cornernodes[i*numberofnodes+j] = curnode;
                else
                    cornernodes[i*numberofnodes+j] = noderenumbering[curnode];
            }
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
    
        reordering[i] = gentools::getreordering(renumbering[i]);    
        // Renumber and reorder the elements in all containers:
        renumber(i, renumbering[i]);
        reorder(i, reordering[i]);  
    }
    
    targetpt->getindisjointregions(indisjointregion);
}

void elements::merge(elements* elstomerge, std::vector<std::vector<int>>& renumbering, std::vector<int>& numduplicates)
{
    std::vector<int> numineachtype = count(); // curvature nodes are included
    
    // Merge the nodes:
    int numnodestomerge = elstomerge->count(0);
    mynodes->setnumber(numineachtype[0] + numnodestomerge - numduplicates[0]);
    
    double* ncs = mynodes->getcoordinates()->data();
    double* ncstomerge = elstomerge->mynodes->getcoordinates()->data();

    for (int i = 0; i < numnodestomerge; i++)
    {
        int renum = renumbering[0][i];
        
        // If the node is already there:
        if (renum < numineachtype[0])
            continue;
            
        ncs[3*renum+0] = ncstomerge[3*i+0];
        ncs[3*renum+1] = ncstomerge[3*i+1];
        ncs[3*renum+2] = ncstomerge[3*i+2];
    }
    
    // Merge the other elements:
    for (int i = 1; i < 8; i++)
    {
        int ne = elstomerge->count(i);
        int nd = numduplicates[i];
    
        // Subelement numbers are below 4:
        for (int j = 0; j < 4; j++)
        {
            int ns = numberofsubelementsineveryelement[i][j];
            
            if (i == j || ns == 0)
                continue;
        
            subelementsinelements[i][j].resize(subelementsinelements[i][j].size() + (ne-nd)*ns);
         
            for (int k = 0; k < ne; k++)
            {
                int renum = renumbering[i][k];
                
                // If the element is already there:
                if (renum < numineachtype[i])
                    continue;
                    
                for (int s = 0; s < ns; s++)
                    subelementsinelements[i][j][renum*ns+s] = renumbering[j][elstomerge->subelementsinelements[i][j][k*ns+s]];
            }
        }
    }
    
    // These containers have been invalidated by adding new elements to this object:
    cleancoordinatedependentcontainers();
    adressedgesatnodes = {};
    edgesatnodes = {};
    adresscellsattype = std::vector<std::vector<int>>(4, std::vector<int>(0));
    cellsattype = std::vector<std::vector<int>>(4, std::vector<int>(0));
}

void elements::merge(std::vector<int> intersectionphysregs, elements* elstomerge)
{
    int meshdim = getdimension();
    std::vector<std::vector<std::vector<int>>*> elementlists(intersectionphysregs.size());
    for (int i = 0; i < intersectionphysregs.size(); i++)
        elementlists[i] = myphysicalregions->get(intersectionphysregs[i])->getelementlist();

    std::vector<int> numineachtype = count(); // curvature nodes are included
    std::vector<int> numduplicates(8, 0);
    // Renumbering of each element to merge:
    std::vector<std::vector<int>> renumberings(8, std::vector<int>(0));
    
    for (int i = 0; i < 8; i++)
    {   
        int numelstomerge = elstomerge->count(i);
        if (numelstomerge == 0)
            continue;
     
        element el(i);
        int eldim = el.getelementdimension();
        // No duplicate check for cells:
        if (eldim >= meshdim)
        {
            renumberings[i] = gentools::getequallyspaced(numineachtype[i], 1, numelstomerge);
            continue;
        }

        // Get the barycenters in this object for the current element type:
        std::vector<bool> isinelementlist;
        int cnt = istypeinelementlists(i, elementlists, isinelementlist, true);
        std::vector<int> targetelemnums;
        gentools::find(isinelementlist, cnt, targetelemnums);
        
        int numtargetels = targetelemnums.size();
        
        // Bring the target barys and the barys to merge in one vector:
        std::vector<double> allcoords(3*(numtargetels+numelstomerge));

        getbarycenters(i, targetelemnums, allcoords.data());
        std::vector<int> allelems = gentools::getequallyspaced(0, 1, numelstomerge);
        elstomerge->getbarycenters(i, allelems, &allcoords[3*numtargetels]);
        
        // Create the renumbering of the elements to merge:
        std::vector<int> rv;
        int numnondupltomerge = gentools::removeduplicates(allcoords, rv) - numtargetels;
        numduplicates[i] = numelstomerge - numnondupltomerge;
        std::vector<int> assignednum(numtargetels+numelstomerge, -1);
        for (int j = 0; j < numtargetels; j++)
            assignednum[rv[j]] = targetelemnums[j];
        int index = 0;
        for (int j = 0; j < numelstomerge; j++)
        {
            int renum = rv[numtargetels+j];
            if (assignednum[renum] == -1)
            {
                assignednum[renum] = numineachtype[i]+index;
                index++;
            }
        }
        
        renumberings[i] = std::vector<int>(numelstomerge);
        for (int j = 0; j < numelstomerge; j++)
            renumberings[i][j] = assignednum[rv[numtargetels+j]];
    }
    
    merge(elstomerge, renumberings, numduplicates);
}

