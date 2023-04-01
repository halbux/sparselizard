#include "dtracker.h"
#include "nodes.h"
#include "elements.h"
#include "physicalregions.h"


bool dtracker::isdefined(void)
{
    return (mynumoverlaplayers >= 0);
}

void dtracker::errorundefined(void)
{
    if (mynumoverlaplayers < 0)
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: DDM context has not been provided" << std::endl;
        log.error();
    }
}

dtracker::dtracker(std::shared_ptr<rawmesh> rm, int globalgeometryskin, int numoverlaplayers)
{
    if (rm == NULL)
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide a NULL rawmesh pointer" << std::endl;
        log.error();
    }
    
    physicalregions* prs = rm->getphysicalregions();
    int meshdim = rm->getmeshdimension();
    
    // Set globalgeometryskin to -1 if empty:
    if (prs->getindex(globalgeometryskin) == -1 || prs->get(globalgeometryskin)->countelements() == 0)
        globalgeometryskin = -1;
    
    if (globalgeometryskin >= 0 && prs->get(globalgeometryskin)->getelementdimension() != meshdim-1)
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: expected " << meshdim-1 << "D elements in the global geometry skin region but found " << prs->get(globalgeometryskin)->getelementdimension() << "D elements (use -1 if empty)" << std::endl;
        log.error();
    }

    myrawmesh = rm;

    myglobalgeometryskin = globalgeometryskin;
    mynumoverlaplayers = numoverlaplayers;
}

std::shared_ptr<rawmesh> dtracker::getrawmesh(void)
{
    if (myrawmesh.expired() == false)
        return myrawmesh.lock();
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot get rawmesh object (weak pointer is expired)" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

void dtracker::setrawmesh(std::shared_ptr<rawmesh> rm)
{
    myrawmesh = rm;
}

bool dtracker::isoverlap(void)
{
    if (mynumoverlaplayers >= 0)
        return (mynumoverlaplayers > 0);
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: number of overlap layers has not been defined (use 0 for no-overlap)" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> dtracker::discoversomeneighbours(int numtrialelements, std::vector<double>& interfaceelembarys, std::vector<int>& neighboursfound)
{ 
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    neighboursfound = {};
    
    int numelementsininterface = interfaceelembarys.size()/3;
    
    std::vector<double> trialbarys(3*numtrialelements, 0.0); // must have this size in any case, barycenters for not alive ranks will be removed below
    if (numelementsininterface > 0)
        gentools::pickcandidates(numtrialelements, interfaceelembarys, trialbarys); // ok if not unique
    
    // Push this in the mpi call:
    trialbarys.push_back(gentools::exactinttodouble(std::min(numelementsininterface,1)));
    
    std::vector<double> alltrialbarys;
    slmpi::allgather(trialbarys, alltrialbarys);
    
    std::vector<double> allisalive = gentools::extract(alltrialbarys, 3*numtrialelements+1, 3*numtrialelements);
    
    // Remove the barycenters of this rank and ranks that are no alive anymore (to avoid removing duplicates with arbitrary all 0 coordinates):
    int numalive = 0;
    for (int r = 0; r < numranks; r++)
    {
        if (r != rank && allisalive[r] == 1.0)
        {
            for (int i = 0; i < numtrialelements; i++)
            {
                alltrialbarys[3*numalive*numtrialelements+3*i+0] = alltrialbarys[3*r*numtrialelements+3*i+0];
                alltrialbarys[3*numalive*numtrialelements+3*i+1] = alltrialbarys[3*r*numtrialelements+3*i+1];
                alltrialbarys[3*numalive*numtrialelements+3*i+2] = alltrialbarys[3*r*numtrialelements+3*i+2];
            }
            numalive++;
        }
    }
    alltrialbarys.resize(3*numalive*numtrialelements);
    numalive += (int)allisalive[rank];
    
    // Return empty if no rank has elements in the interface:
    if (numalive == 0)
        return {};

    // Find trial elements coordinates on the interface:  
    std::vector<int> isfound;
    if (numelementsininterface > 0)
        gentools::findcoordinates(interfaceelembarys, alltrialbarys, isfound);
        
    // Deduce which are neighbours. To avoid additional mpi length exchange steps the number of neighbours is limited.
    std::vector<int> neighbourlist(numtrialelements, -1);
    if (numelementsininterface > 0)
    {
        int index = 0, rindex = 0; // rindex to take into account the ranks with removed barycenters
        std::vector<bool> isalreadythere(numranks, false);
        for (int r = 0; r < numranks; r++)
        {
            if (r != rank && allisalive[r] == 1.0)
            {
                for (int i = 0; i < numtrialelements; i++)
                {
                    if (index < numtrialelements && isalreadythere[r] == false && isfound[rindex*numtrialelements+i] != -1)
                    {
                        isalreadythere[r] = true;
                        neighbourlist[index] = r;
                        index++;
                    }
                }

                rindex++;
            }
        }
    }
    
    // Push this in the mpi call:
    neighbourlist.push_back(numelementsininterface);
    
    std::vector<int> allneighbourlist;
    slmpi::allgather(neighbourlist, allneighbourlist);
    
    std::vector<int> allnumelementsininterface = gentools::extract(allneighbourlist, numtrialelements+1, numtrialelements);
    
    // Deduce the neighbours for the current rank:
    if (numelementsininterface > 0)
    {
        std::vector<bool> isneighbour(numranks, false);
        for (int r = 0; r < numranks; r++)
        {
            for (int i = 0; i < numtrialelements; i++)
            {
                int cur = allneighbourlist[r*numtrialelements + i];
                if (cur < 0)
                    continue;
            
                // Neighbour relations are symmetric:
                if (rank == r)
                    isneighbour[cur] = true;
                else
                {
                    if (rank == cur)
                        isneighbour[r] = true;
                }
            }
        }
    
        for (int r = 0; r < numranks; r++)
        {
             if (isneighbour[r])
                neighboursfound.push_back(r); // push_back here is ok (low number of neighbours)
        }
    }
    
    return allnumelementsininterface;
}

void dtracker::discoverinterfaces(std::vector<int> neighbours, std::vector<double>& interfaceelembarys, std::vector<int>& allnumelementsininterface, std::vector<int>& inneighbour)
{
    int rank = slmpi::getrank();
    
    int numelementsininterface = interfaceelembarys.size()/3;
    
    inneighbour = std::vector<int>(numelementsininterface, -1);
    
    int numneighbours = neighbours.size();
    if (numneighbours == 0)
        return;
        
    // Split neighbours into lower ranks and higher ranks:
    std::vector<int> lowerranks, higherranks;
    gentools::split(neighbours, rank, lowerranks, higherranks);
    
    // Send all interface barycenter coordinates to the higher ranks and receive from the lower ranks in a vector:
    int totnumcand = 0;
    for (int n = 0; n < lowerranks.size(); n++)
        totnumcand += allnumelementsininterface[lowerranks[n]];
        
    std::vector<double> candidatebarys(3*totnumcand);
    
    // It is safer not to send with the same buffer:
    std::vector<std::vector<double>> sends(higherranks.size(), interfaceelembarys);
    
    std::vector<int> sendlens(numneighbours, 0);
    std::vector<double*> sendbuffers(numneighbours, NULL);
    std::vector<int> receivelens(numneighbours, 0);
    std::vector<double*> receivebuffers(numneighbours, NULL);
    
    int pos = 0;
    if (totnumcand > 0)
    {
        for (int n = 0; n < lowerranks.size(); n++)
        {
            int len = 3*allnumelementsininterface[lowerranks[n]];
            receivelens[n] = len;
            receivebuffers[n] = &candidatebarys[pos];
            pos += len;
        }
    }
    for (int n = 0; n < higherranks.size(); n++)
    {
        sendlens[lowerranks.size()+n] = 3*numelementsininterface;
        sendbuffers[lowerranks.size()+n] = sends[n].data();
    }
    
    slmpi::exchange(neighbours, sendlens, sendbuffers, receivelens, receivebuffers);
    
    // Find matches:
    std::vector<int> posfound;
    gentools::findcoordinates(interfaceelembarys, candidatebarys, posfound);
    
    std::vector<std::vector<bool>> isfound(lowerranks.size());
    
    pos = 0;
    for (int n = 0; n < lowerranks.size(); n++)
    {
        int len = allnumelementsininterface[lowerranks[n]];
        
        isfound[n] = std::vector<bool>(len, false);
        
        for (int i = 0; i < len; i++)
        {
            int pf = posfound[pos+i];
            if (pf != -1)
            {
                inneighbour[pf] = lowerranks[n];
                isfound[n][i] = true;
            }
        }
        pos += len;
    }
    
    // Send match information to every neighbour of lower rank:
    std::vector<std::vector<int>> tosend(numneighbours, std::vector<int>(0));
    for (int n = 0; n < lowerranks.size(); n++)
        gentools::pack(isfound[n], tosend[n]);

    std::vector<std::vector<int>> toreceive(numneighbours, std::vector<int>(0));
    for (int n = 0; n < higherranks.size(); n++)
        toreceive[lowerranks.size()+n].resize(gentools::getpackedsize(numelementsininterface));

    slmpi::exchange(neighbours, tosend, toreceive);

    for (int n = lowerranks.size(); n < numneighbours; n++)
    {
        std::vector<bool> wasfound;
        gentools::unpack(numelementsininterface, toreceive[n], wasfound);
    
        for (int i = 0; i < numelementsininterface; i++)
        {
            if (wasfound[i])
                inneighbour[i] = neighbours[n];
        }
    }
}

bool dtracker::discovercrossinterfaces(std::vector<int>& interfacenodelist, std::vector<int>& interfaceedgelist, std::vector<std::vector<bool>>& isnodeinneighbours, std::vector<std::vector<bool>>& isedgeinneighbours)
{
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();

    int numranks = slmpi::count();
    
    int numnodes = nds->count();
    int numedges = els->count(1);
    
    std::vector<int> neighbours = {};
    for (int i = 0; i < numranks; i++)
    {
        if (isnodeinneighbours[i].size() > 0 || isedgeinneighbours[i].size() > 0)
        {
            // Make sure both are defined together:
            if (isnodeinneighbours[i].size() == 0)
                isnodeinneighbours[i] = std::vector<bool>(numnodes, false);
            if (isedgeinneighbours[i].size() == 0)
                isedgeinneighbours[i] = std::vector<bool>(numedges, false);
        
            neighbours.push_back(i);
        }
    }
    int numneighbours = neighbours.size();
    
    std::vector<double>* ncs = nds->getcoordinates();
    std::vector<double>* edgebarys = els->getbarycenters(1);
    
    // For each neighbour pair make a container that stores all nodes/edges in the interface intersection.
    // Neighbour relations are symmetric thus only indexes i*numneighbours+j with j >= i+1 are considered.
    std::vector<int> numnodesinneighbourpair(numneighbours*numneighbours, 0);
    std::vector<int> numedgesinneighbourpair(numneighbours*numneighbours, 0);
    
    std::vector<std::vector<bool>> nodesinneighbourpair(numneighbours*numneighbours, std::vector<bool>(0));
    std::vector<std::vector<bool>> edgesinneighbourpair(numneighbours*numneighbours, std::vector<bool>(0));
    
    // Treat neighbour pair (m,n):
    for (int n = 0; n < numneighbours; n++)
    {
        // Start at n+1 because neighbour relations are symmetric (counted twice otherwise) and not to itself:
        for (int m = n+1; m < numneighbours; m++)
        {
            for (int i = 0; i < numnodes; i++)
            {
                if (isnodeinneighbours[neighbours[n]][i] && isnodeinneighbours[neighbours[m]][i])
                {
                    if (nodesinneighbourpair[n*numneighbours+m].size() == 0)
                        nodesinneighbourpair[n*numneighbours+m] = std::vector<bool>(numnodes, false);
                
                    numnodesinneighbourpair[n*numneighbours+m]++;
                    nodesinneighbourpair[n*numneighbours+m][i] = true;
                }
            }
            for (int i = 0; i < numedges; i++)
            {
                if (isedgeinneighbours[neighbours[n]][i] && isedgeinneighbours[neighbours[m]][i])
                {
                    if (edgesinneighbourpair[n*numneighbours+m].size() == 0)
                        edgesinneighbourpair[n*numneighbours+m] = std::vector<bool>(numedges, false);
                
                    numedgesinneighbourpair[n*numneighbours+m]++;
                    edgesinneighbourpair[n*numneighbours+m][i] = true;
                }
            }
        }
    }
    
    // Package node and edge barycenters for both neighbours in each neighbour pair:
    std::vector<std::vector<double>> packaged(numneighbours*numneighbours, std::vector<double>(0));
    for (int i = 0; i < numneighbours*numneighbours; i++)
    {
        int totnumtosend = numnodesinneighbourpair[i] + numedgesinneighbourpair[i];
        
        if (totnumtosend > 0)
        {
            packaged[i] = std::vector<double>(3*totnumtosend+2);
            packaged[i][0] = 3*gentools::exactinttodouble(numnodesinneighbourpair[i]);
            packaged[i][1] = 3*gentools::exactinttodouble(numedgesinneighbourpair[i]);
            
            gentools::selectcoordinates(nodesinneighbourpair[i], *ncs, &(packaged[i][2]));
            gentools::selectcoordinates(edgesinneighbourpair[i], *edgebarys, &(packaged[i][2+3*numnodesinneighbourpair[i]]));
        }
    }
    
    std::vector<std::vector<double>> dataforeachneighbour;
    gentools::pack(neighbours, packaged, dataforeachneighbour);
    
    // Exchange all interface barycenter coordinates with every neighbour:
    std::vector<int> sendlens(numneighbours);
    std::vector<double*> sendbuffers(numneighbours);
    std::vector<int> receivelens(numneighbours);
    std::vector<double*> receivebuffers(numneighbours);
    
    std::vector<std::vector<double>> datafromeachneighbour(numneighbours, std::vector<double>(0));
    
    for (int n = 0; n < numneighbours; n++)
    {
        int datasize = dataforeachneighbour[n].size();
        sendlens[n] = datasize;
        sendbuffers[n] = dataforeachneighbour[n].data();
    }
    slmpi::exchange(neighbours, sendlens, receivelens);
    
    for (int n = 0; n < numneighbours; n++)
    {
        datafromeachneighbour[n].resize(receivelens[n]);
        receivebuffers[n] = datafromeachneighbour[n].data();
    }
    slmpi::exchange(neighbours, sendlens, sendbuffers, receivelens, receivebuffers);
    
    // Unpack:
    std::vector< std::vector<std::vector<double>> > unpackedcoords(numneighbours);
    std::vector<std::vector<int>> candidateneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
        candidateneighbours[n] = gentools::unpack(datafromeachneighbour[n], unpackedcoords[n]);
    
    std::vector<double> allreceivednodecoords, allreceivededgecoords;
    std::vector<std::vector<int>> lennodedataingroup, lenedgedataingroup;
    gentools::split(unpackedcoords, allreceivednodecoords, allreceivededgecoords, lennodedataingroup, lenedgedataingroup);
    
    std::vector<double> interfacenodesbarys;
    els->getbarycenters(0, interfacenodelist, interfacenodesbarys);
    std::vector<double> interfaceedgesbarys;
    els->getbarycenters(1, interfaceedgelist, interfaceedgesbarys);
    
    std::vector<int> posnodefound;
    int nnf = gentools::findcoordinates(interfacenodesbarys, allreceivednodecoords, posnodefound);
    std::vector<int> posedgefound;
    int nef = gentools::findcoordinates(interfaceedgesbarys, allreceivededgecoords, posedgefound);

    // All nodes and edges must have been found or something went wrong:
    if (nnf != posnodefound.size() || nef != posedgefound.size())
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: algorithm to detect cross-interfaces did not succeed (at least one node or edge barycenter was not matched on a target neighbour)" << std::endl;
        log.error();
    }
    
    // Update 'isnodeinneighbours' and 'isedgeinneighbours':
    int nodeindex = 0;
    int edgeindex = 0;
    
    int isanynewadded = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < unpackedcoords[n].size(); i++)
        {
            int nn = lennodedataingroup[n][i]/3;
            int ne = lenedgedataingroup[n][i]/3;
            
            if (nn == 0 && ne == 0)
                continue;
            
            int curcanneighour = candidateneighbours[n][i];
            
            if (isnodeinneighbours[curcanneighour].size() == 0)
                isnodeinneighbours[curcanneighour] = std::vector<bool>(numnodes, false);
            if (isedgeinneighbours[curcanneighour].size() == 0)
                isedgeinneighbours[curcanneighour] = std::vector<bool>(numedges, false);

            for (int k = 0; k < nn; k++)
            {
                int nodetoadd = interfacenodelist[posnodefound[nodeindex+k]];
                if (isnodeinneighbours[curcanneighour][nodetoadd] == false)
                {
                    isnodeinneighbours[curcanneighour][nodetoadd] = true;
                    isanynewadded = 1;
                }
            }
            for (int k = 0; k < ne; k++)
            {
                int edgetoadd = interfaceedgelist[posedgefound[edgeindex+k]];
      
                if (isedgeinneighbours[curcanneighour][edgetoadd] == false)
                {
                    isedgeinneighbours[curcanneighour][edgetoadd] = true;
                    isanynewadded = 1;
                }
            }
        
            nodeindex += nn;
            edgeindex += ne;
        }
    }

    // Status must be checked on all ranks to decide:
    slmpi::sum(1, &isanynewadded);
    
    return (isanynewadded != 0);
}

void dtracker::defineinneroverlaps(void)
{
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int numranks = slmpi::count();

    myinneroverlaps = std::vector<int>(numranks, -1);
    myinneroverlapskins = std::vector<int>(numranks, -1);

    for (int n = 0; n < myneighbours.size(); n++)
    {
        regiondefiner regdef(*nds, *els, *prs);

        int cn = myneighbours[n];
        int firstnewpr = prs->getmaxphysicalregionnumber()+1;

        myinneroverlaps[cn] = firstnewpr;
        myinneroverlapskins[cn] = firstnewpr+1;

        for (int dim = 0; dim < 3; dim++)
        {
            int ci = mynooverlapinterfaces[3*cn+dim];
            if (ci >= 0)
                regdef.regionlayer(myinneroverlaps[cn], -1, ci, mynumoverlaplayers);
        }

        regdef.regionskin(myinneroverlapskins[cn], myinneroverlaps[cn]);

        // Duplicated cells in the inner overlap are automatically removed:
        regdef.defineregions();
    }
}

void dtracker::exchangeoverlaps(void)
{
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int numranks = slmpi::count();
    
    double* ncs = nds->getcoordinates()->data();
    int curvatureorder = els->getcurvatureorder();

    int numneighbours = myneighbours.size();

    std::vector<int> numelemsnooverlap = els->count();
    
    // Node numbers and node coordinates for all cells in the inner overlaps:
    std::vector<std::vector<int>> nodenumsforeachneighbour(numneighbours);
    std::vector<std::vector<double>> coordsforeachneighbour(numneighbours);
        
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<std::vector<int>>* inneroverlapcells = prs->get(myinneroverlaps[cn])->getelementlist();
        
        int prealloc = 0;
        for (int i = 0; i < 8; i++)
        {
            element el(i, curvatureorder);
            prealloc += inneroverlapcells->at(i).size() * el.countcurvednodes();
        }
        nodenumsforeachneighbour[n] = std::vector<int>(prealloc);
        
        // Get the nodes in each cell:
        int ni = 0;
        for (int i = 0; i < 8; i++)
        {
            element el(i, curvatureorder);
            int ncn = el.countcurvednodes();
            
            for (int j = 0; j < inneroverlapcells->at(i).size(); j++)
            {
                int elem = inneroverlapcells->at(i)[j];
                for (int k = 0; k < ncn; k++)
                    nodenumsforeachneighbour[n][ni+k] = els->getsubelement(0, i, elem, k);
                ni += ncn;
            }
        }
        
        // Get the renumbering to make the node numbers consecutive:
        std::vector<int> renumvec;
        int numuniques = gentools::squeeze(nodenumsforeachneighbour[n], nds->count(), renumvec);
        
        // Make the nodes consecutive and gather all coordinates to send:
        coordsforeachneighbour[n] = std::vector<double>(3*numuniques);
        for (int i = 0; i < prealloc; i++)
        {
            int curnode = nodenumsforeachneighbour[n][i];
            int renumbered = renumvec[curnode];
            
            nodenumsforeachneighbour[n][i] = renumbered;
            coordsforeachneighbour[n][3*renumbered+0] = ncs[3*curnode+0];
            coordsforeachneighbour[n][3*renumbered+1] = ncs[3*curnode+1];
            coordsforeachneighbour[n][3*renumbered+2] = ncs[3*curnode+2];
        }
    }
            
    // Exchange the node numbers and coordinates with each neighbour:
    std::vector<int> sendlens(11*numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        sendlens[11*n+0] = nodenumsforeachneighbour[n].size();
        sendlens[11*n+1] = coordsforeachneighbour[n].size();
        // Send the number of inner overlap cells of each type as well as the curvature order:
        std::vector<std::vector<int>>* inneroverlapcells = prs->get(myinneroverlaps[myneighbours[n]])->getelementlist();
        for (int i = 0; i < 8; i++)
            sendlens[11*n+2+i] = inneroverlapcells->at(i).size();
        sendlens[11*n+10] = curvatureorder;
    }
    std::vector<int> reclens;
    slmpi::exchange(myneighbours, sendlens, reclens);

    std::vector<std::vector<int>> nodenumsfromeachneighbour(numneighbours);
    std::vector<std::vector<double>> coordsfromeachneighbour(numneighbours);
    
    for (int n = 0; n < numneighbours; n++)
    {
        nodenumsfromeachneighbour[n].resize(reclens[11*n+0]);
        coordsfromeachneighbour[n].resize(reclens[11*n+1]);
    }

    slmpi::exchange(myneighbours, nodenumsforeachneighbour, nodenumsfromeachneighbour);
    slmpi::exchange(myneighbours, coordsforeachneighbour, coordsfromeachneighbour); 

    // Make sure the curvature order is identical on all neighbours:
    for (int n = 0; n < numneighbours; n++)
    {
        if (reclens[11*n+10] != curvatureorder)
        {
            logs log;
            log.msg() << "Error in 'dtracker' object: expected the same curvature order on all mesh parts" << std::endl;
            log.error();
        }
    }

    // Populate the temporary 'nodes' object:
    nodes tmpnds;
    
    int totnumrecnodes = 0;
    for (int n = 0; n < numneighbours; n++)
        totnumrecnodes += coordsfromeachneighbour[n].size()/3;
    
    tmpnds.setnumber(totnumrecnodes);
    double* tmpncs = tmpnds.getcoordinates()->data();
    
    int cindex = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < coordsfromeachneighbour[n].size(); i++)
            tmpncs[cindex+i] = coordsfromeachneighbour[n][i];
        cindex += coordsfromeachneighbour[n].size();
    }
    
    // Populate the temporary 'elements' object:
    disjointregions tmpdrs;
    physicalregions tmpprs(tmpdrs);
    elements tmpels(tmpnds, tmpprs, tmpdrs);
    
    int nodeshift = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int ni = 0;
        for (int i = 0; i < 8; i++)
        {
            element el(i, curvatureorder);
            int ncn = el.countcurvednodes();
            std::vector<int> nodelist(ncn);
            
            int numelemsintype = reclens[11*n+2+i];
            for (int j = 0; j < numelemsintype; j++)
            {
                for (int k = 0; k < ncn; k++)
                    nodelist[k] = nodenumsfromeachneighbour[n][ni+k] + nodeshift;
                ni += ncn;
                
                tmpels.add(i, curvatureorder, nodelist);
            }
        }
        nodeshift += coordsfromeachneighbour[n].size()/3;
    }
    if (numneighbours > 0)
        tmpels.explode();
    
    // Merge while removing the duplicates:
    std::vector<int> nooverlapregs = {};
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        for (int dim = 0; dim < 3; dim++)
        {
            int cr = mynooverlapinterfaces[3*cn+dim];
            if (cr >= 0)
                nooverlapregs.push_back(cr);
        }
    }
    if (numneighbours > 0)
        els->merge(nooverlapregs, &tmpels);
    
    // Define the outer overlap regions and their skins:
    myouteroverlaps = std::vector<int>(numranks, -1);
    myouteroverlapskins = std::vector<int>(numranks, -1);

    std::vector<int> offsets = numelemsnooverlap;
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];

        int firstnewpr = prs->getmaxphysicalregionnumber()+1;
        myouteroverlaps[cn] = firstnewpr;
        myouteroverlapskins[cn] = firstnewpr+1;

        physicalregion* oopr = prs->get(firstnewpr);

        for (int i = 0; i < 8; i++)
        {
            int numelemsintype = reclens[11*n+2+i];

            for (int j = 0; j < numelemsintype; j++)
                oopr->addelement(i, offsets[i]+j); // outer overlap cells have been appended to 'els'

            offsets[i] += numelemsintype;
        }

        regiondefiner regdef(*nds, *els, *prs);
        regdef.regionskin(myouteroverlapskins[cn], myouteroverlaps[cn]);
        regdef.defineregions();
    }
}

void dtracker::exchangephysicalregions(void)
{
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int rank = slmpi::getrank();
    
    int numneighbours = myneighbours.size();

    std::vector<int> ddmregs = listddmregions();
    std::vector<bool> isddmpr(prs->getmaxphysicalregionnumber()+1, false);
    for (int i = 0; i < ddmregs.size(); i++)
        isddmpr[ddmregs[i]] = true;

    // Get the info of all physical regions in which each element is:
    std::vector<std::vector<int>> addressesinphysreglist(8);
    std::vector<std::vector<int>> physreglist(8);
    for (int i = 0; i < 8; i++)
        prs->inphysicalregions(i, els->count(i), addressesinphysreglist[i], physreglist[i]);

    // Get the list of (sub)elements in each inner and outer overlap region:
    std::vector<std::vector<std::vector<int>>> elemsininneroverlaps(numneighbours, std::vector<std::vector<int>>(8, std::vector<int>(0)));
    std::vector<std::vector<std::vector<int>>> elemsinouteroverlaps(numneighbours, std::vector<std::vector<int>>(8, std::vector<int>(0)));
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<std::vector<int>>* inneroverlapcells = prs->get(myinneroverlaps[cn])->getelementlist();
        std::vector<std::vector<int>>* outeroverlapcells = prs->get(myouteroverlaps[cn])->getelementlist();
        
        for (int i = 0; i < 8; i++)
        {
            els->follow(inneroverlapcells, i, elemsininneroverlaps[n][i]);
            els->follow(outeroverlapcells, i, elemsinouteroverlaps[n][i]);
        }
    }

    std::vector<std::vector<int>> physreglistsforeachneighbour(numneighbours);
    std::vector<int> sendlens(numneighbours), reclens(numneighbours);

    for (int n = 0; n < numneighbours; n++)
    {
        int prealloc = 0;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < elemsininneroverlaps[n][i].size(); j++)
            {
                int curelem = elemsininneroverlaps[n][i][j];
                prealloc += 1 + addressesinphysreglist[i][curelem+1]-addressesinphysreglist[i][curelem];
            }
        }
        
        physreglistsforeachneighbour[n].resize(prealloc);
        
        int pi = 0;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < elemsininneroverlaps[n][i].size(); j++)
            {
                int curelem = elemsininneroverlaps[n][i][j];
                int numphysregsincurelem = addressesinphysreglist[i][curelem+1]-addressesinphysreglist[i][curelem];
                
                int actualnumphysregsincurelem = 0;
                for (int l = 0; l < numphysregsincurelem; l++)
                {
                    int cpr = physreglist[i][addressesinphysreglist[i][curelem]+l];
                    if (isddmpr[cpr] == false)
                    {
                        physreglistsforeachneighbour[n][pi+1+actualnumphysregsincurelem] = cpr;
                        actualnumphysregsincurelem++;
                    }
                }
                physreglistsforeachneighbour[n][pi] = actualnumphysregsincurelem;
                
                pi += 1+actualnumphysregsincurelem;
            }
        }
        physreglistsforeachneighbour[n].resize(pi);
        
        gentools::compresszeros(physreglistsforeachneighbour[n]);
        
        sendlens[n] = physreglistsforeachneighbour[n].size();
    }
    
    slmpi::exchange(myneighbours, sendlens, reclens);
    
    std::vector<std::vector<int>> physreglistsfromeachneighbour(numneighbours);
    for (int n = 0; n < numneighbours; n++)
        physreglistsfromeachneighbour[n].resize(reclens[n]);

    slmpi::exchange(myneighbours, physreglistsforeachneighbour, physreglistsfromeachneighbour);
    
    // Add the physical regions to the elements (the outer overlap matches the inner overlap for each neighbour):
    std::vector<physicalregion*> allprs(prs->getmaxphysicalregionnumber()+1, NULL);
    for (int n = 0; n < numneighbours; n++)
    {
        gentools::decompresszeros(physreglistsfromeachneighbour[n]);
        
        int pi = 0;
        for (int i = 0; i < 8; i++)
        {
            for (int j = 0; j < elemsinouteroverlaps[n][i].size(); j++)
            {
                int numphysregsincurelem = physreglistsfromeachneighbour[n][pi];
                
                for (int l = 0; l < numphysregsincurelem; l++)
                {
                    int curpreg = physreglistsfromeachneighbour[n][pi+1+l];

                    // The received physical region number might be larger than any in this rank:
                    if (curpreg >= allprs.size())
                    {
                        int prevlen = allprs.size();
                        allprs.resize(curpreg+1);
                        isddmpr.resize(curpreg+1);
                        for (int p = prevlen; p < allprs.size(); p++)
                        {
                            allprs[p] = NULL;
                            isddmpr[p] = false;
                        }
                    }
                    
                    if (allprs[curpreg] == NULL)
                    {
                        if (isddmpr[curpreg])
                        {
                            logs log;
                            log.msg() << "Error in 'dtracker' object: cannot merge physical region " << curpreg << " from the inner overlap of rank " << myneighbours[n] << " into the outer overlap of rank " << rank << " (it is a DDM owned region on the latter rank)" << std::endl;
                            log.error();
                        }
                        
                        allprs[curpreg] = prs->get(curpreg);
                    }
                        
                    allprs[curpreg]->addelement(i, elemsinouteroverlaps[n][i][j]);
                }
                pi += 1+numphysregsincurelem;
            }
        }
    }
    
    // Remove the duplicated elements that have been created above:
    for (int i = 0; i < allprs.size(); i++)
    {
        if (allprs[i] != NULL)
            allprs[i]->removeduplicatedelements();
    }
}

void dtracker::defineouteroverlapinterfaces(void)
{
    // The outer overlap interface with a neighbour is the intersection of:
    //
    // - the skin of the overlapped domain minus the global geometry skin
    // - the outer overlap skin

    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int numranks = slmpi::count();

    int numneighbours = myneighbours.size();

    myouteroverlapinterfaces = std::vector<int>(3*numranks, -1);
    
    if (numneighbours == 0)
        return;
    
    int skinregion = prs->getmaxphysicalregionnumber()+1;
    int candidatesregion = skinregion+1;
    
    regiondefiner regdef(*nds, *els, *prs);

    regdef.regionskin(skinregion, -1);
    // This domain might not be in contact with the global geometry skin:
    if (myglobalgeometryskin >= 0)
        regdef.regionexclusion(candidatesregion, skinregion, {myglobalgeometryskin});
    else
        candidatesregion = skinregion;

    regdef.defineregions();

    std::vector<std::vector<int>>* candidateelementlist = prs->get(candidatesregion)->getelementlist();
    
    int numinterfacetype = els->countinterfaceelementtypes();
    std::vector<std::vector<bool>> iscandidate(numinterfacetype, std::vector<bool>(0));
    for (int i = 0; i < numinterfacetype; i++)
        els->istypeinelementlists(i, {candidateelementlist}, iscandidate[i], false);
    
    // Remove the construction regions:
    prs->remove({skinregion, candidatesregion}, false);
    
    // Create the outer overlap interfaces as the intersection between the above created region and each outer overlap skin:
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<std::vector<bool>> isinoos(numinterfacetype, std::vector<bool>(0));
        for (int i = 0; i < numinterfacetype; i++)
            els->istypeinelementlists(i, {prs->get(myouteroverlapskins[cn])->getelementlist()}, isinoos[i], false);
    
        // Add the intersection faces, edges and nodes (if any):
        std::vector<physicalregion*> dimprs(3, NULL);
        for (int i = numinterfacetype-1; i >= 0; i--)
        {
            element el(i);
            int dim = el.getelementdimension();
            int nn = el.countnodes();
            int ne = el.countedges();
            
            for (int j = 0; j < els->count(i); j++)
            {
                if (iscandidate[i][j] && isinoos[i][j])
                {
                    if (dimprs[dim] == NULL)
                    {
                        myouteroverlapinterfaces[3*cn+dim] = prs->getmaxphysicalregionnumber()+1;
                        dimprs[dim] = prs->get(myouteroverlapinterfaces[3*cn+dim]);
                    }
                    dimprs[dim]->addelement(i, j);
                    
                    // Remove the underlying edges and nodes:
                    if (dim > 1)
                    {
                        for (int k = 0; k < ne; k++)
                            isinoos[1][els->getsubelement(1,i,j,k)] = false;
                    }
                    if (dim > 0)
                    {
                        for (int k = 0; k < nn; k++)
                            isinoos[0][els->getsubelement(0,i,j,k)] = false;
                    }
                }
            }
        }
    }
}

void dtracker::defineinneroverlapinterfaces(void)
{
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int numranks = slmpi::count();

    int numneighbours = myneighbours.size();
    
    myinneroverlapinterfaces = std::vector<int>(3*numranks, -1);
    
    if (numneighbours == 0)
        return;
    
    int numinterfacetype = els->countinterfaceelementtypes();
    
    // Get the list of (sub)elements in each inner and outer overlap skin region:
    std::vector<std::vector<std::vector<int>>> elemsininneroverlapskins(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    std::vector<std::vector<std::vector<int>>> elemsinouteroverlapskins(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<std::vector<int>>* inneroverlapelems = prs->get(myinneroverlaps[cn])->getelementlist();
        std::vector<std::vector<int>>* outeroverlapelems = prs->get(myouteroverlaps[cn])->getelementlist();
        
        std::vector<std::vector<int>>* inneroverlapskinelems = prs->get(myinneroverlapskins[cn])->getelementlist();
        std::vector<std::vector<int>>* outeroverlapskinelems = prs->get(myouteroverlapskins[cn])->getelementlist();
        
        for (int i = 0; i < numinterfacetype; i++)
        {
            els->follow(inneroverlapelems, i, elemsininneroverlapskins[n][i], {inneroverlapskinelems});
            els->follow(outeroverlapelems, i, elemsinouteroverlapskins[n][i], {outeroverlapskinelems});
        }
    }
    
    // Send bool data as ints:
    std::vector<int> iosnumbits(numneighbours, 0);
    std::vector<std::vector<int>> dataforneighbours(numneighbours), datafromneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        int oosnumbits = 0;
        for (int i = 0; i < numinterfacetype; i++)
            oosnumbits += elemsinouteroverlapskins[n][i].size();
            
        std::vector<bool> isininterface(oosnumbits, false);
        
        int pos = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            element el(i);
            int dim = el.getelementdimension();
        
            int ne = elemsinouteroverlapskins[n][i].size();
            int cr = myouteroverlapinterfaces[3*cn+dim];
            
            if (ne > 0 && cr >= 0)
            {
                std::vector<bool> issubinooi;
                els->istypeinelementlists(i, {prs->get(cr)->getelementlist()}, issubinooi, false);
                
                for (int j = 0; j < ne; j++)
                {
                    if (issubinooi[elemsinouteroverlapskins[n][i][j]])
                        isininterface[pos+j] = true;
                }
            }
            pos += ne;
        }
        
        gentools::pack(isininterface, dataforneighbours[n]);
        
        for (int i = 0; i < numinterfacetype; i++)
            iosnumbits[n] += elemsininneroverlapskins[n][i].size();
        datafromneighbours[n].resize(gentools::getpackedsize(iosnumbits[n]));
    }
    
    slmpi::exchange(myneighbours, dataforneighbours, datafromneighbours);
        
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<bool> isininterface;
        gentools::unpack(iosnumbits[n], datafromneighbours[n], isininterface);
    
        std::vector<physicalregion*> dimprs(3, NULL);
        
        int pos = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            element el(i);
            int dim = el.getelementdimension();
            
            int ne = elemsininneroverlapskins[n][i].size();
            for (int j = 0; j < ne; j++)
            {
                if (isininterface[pos+j])
                {
                    if (dimprs[dim] == NULL)
                    {
                        myinneroverlapinterfaces[3*cn+dim] = prs->getmaxphysicalregionnumber()+1;
                        dimprs[dim] = prs->get(myinneroverlapinterfaces[3*cn+dim]);
                    }
                    dimprs[dim]->addelement(i, elemsininneroverlapskins[n][i][j]);
                }
            }
            pos += ne;
        }
    }
}

void dtracker::mapnooverlapinterfaces(void)
{
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    
    int numneighbours = myneighbours.size();
    
    int numinterfacetype = els->countinterfaceelementtypes();
    
    // Create a container listing all subelements in each no-overlap interface:
    std::vector<int> numsubsoneachneighbour(numneighbours, 0);
    std::vector<std::vector<std::vector<int>>> subsininterfaces(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        for (int i = 0; i < numinterfacetype; i++)
        {   
            std::vector<std::vector<std::vector<int>>*> interfaceelems(3, NULL);
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = mynooverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                    interfaceelems[dim] = prs->get(cr)->getelementlist();
            }
            std::vector<bool> isininterface;
            int numininterface = els->istypeinelementlists(i, interfaceelems, isininterface, false); // no curvature nodes
            gentools::find(isininterface, numininterface, subsininterfaces[n][i]);
            
            numsubsoneachneighbour[n] += numininterface;
        }
    }
    
    // Split neighbours into lower ranks and higher ranks (neighbours are sorted ascendingly):
    std::vector<int> lowerranks, higherranks;
    gentools::split(myneighbours, rank, lowerranks, higherranks);
    
    // Prepare all barycenters to send to the higher ranks:
    std::vector<std::vector<double>> barysforeachneighbour(numneighbours, std::vector<double>(0));
    for (int n = lowerranks.size(); n < numneighbours; n++)
    {
        barysforeachneighbour[n] = std::vector<double>(3*numsubsoneachneighbour[n]);
        
        if (numsubsoneachneighbour[n] > 0)
        {
            int pos = 0;
            for (int i = 0; i < numinterfacetype; i++)
            {
                els->getbarycenters(i, subsininterfaces[n][i], &barysforeachneighbour[n][pos]);
                pos += 3*subsininterfaces[n][i].size();
            }
        }
    }
    std::vector<std::vector<double>> barysfromeachneighbour(numneighbours, std::vector<double>(0));
    for (int n = 0; n < lowerranks.size(); n++)
        barysfromeachneighbour[n].resize(3*numsubsoneachneighbour[n]);
    
    slmpi::exchange(myneighbours, barysforeachneighbour, barysfromeachneighbour);
    
    // Treat all barycenters sent by the lower ranks:
    std::vector<std::vector<std::vector<int>>> posfounds(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    for (int n = 0; n < lowerranks.size(); n++)
    {
        int pos = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            int num = subsininterfaces[n][i].size();
            std::vector<double> recbarys(3*num);
            for (int j = 0; j < 3*num; j++)
                recbarys[j] = barysfromeachneighbour[n][pos+j];
            
            pos += 3*num;
            
            std::vector<double> targetbarys;
            els->getbarycenters(i, subsininterfaces[n][i], targetbarys);
            
            int numfound = gentools::findcoordinates(targetbarys, recbarys, posfounds[n][i]);
            if (numfound != num)
            {
                logs log;
                log.msg() << "Error in 'dtracker' object: no-overlap mapping failed because some interface elements could not be matched" << std::endl;
                log.error();
            }
        }
    }
    
    // Exchange the subs numbers for mapping:
    std::vector<std::vector<int>> elemnumbersforeachneighbour(numneighbours), elemnumbersfromeachneighbour(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        // Also include the number of interface elements of each type (needed to preallocate the map):
        elemnumbersforeachneighbour[n].resize(numsubsoneachneighbour[n] + numinterfacetype);
        elemnumbersfromeachneighbour[n].resize(numsubsoneachneighbour[n] + numinterfacetype);
    
        int pos = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            int num = subsininterfaces[n][i].size();
            
            // Lower rank neighbour:
            if (n < lowerranks.size())
            {
                for (int j = 0; j < num; j++)
                    elemnumbersforeachneighbour[n][pos+j] = subsininterfaces[n][i][posfounds[n][i][j]];
            }
            else
            {
                for (int j = 0; j < num; j++)
                    elemnumbersforeachneighbour[n][pos+j] = subsininterfaces[n][i][j];
            }
            pos += num;
        }
        for (int i = 0; i < numinterfacetype; i++)
            elemnumbersforeachneighbour[n][pos+i] = els->count(i);
    }
    
    slmpi::exchange(myneighbours, elemnumbersforeachneighbour, elemnumbersfromeachneighbour);
    
    // Create and populate the map:
    mymaptothisdomain = std::vector<std::vector<std::vector<int>>>(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    
    for (int n = 0; n < numneighbours; n++)
    {
        // Preallocate the map:
        for (int i = 0; i < numinterfacetype; i++)
            mymaptothisdomain[n][i] = std::vector<int>(elemnumbersfromeachneighbour[n][elemnumbersfromeachneighbour[n].size()-numinterfacetype+i], -1);
    
        int pos = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            int num = subsininterfaces[n][i].size();
            
            // Lower rank neighbour:
            if (n < lowerranks.size())
            {
                for (int j = 0; j < num; j++)
                {
                    int recnum = elemnumbersfromeachneighbour[n][pos+j];
                    mymaptothisdomain[n][i][recnum] = subsininterfaces[n][i][posfounds[n][i][j]];
                }
            }
            else
            {
                for (int j = 0; j < num; j++)
                {
                    int recnum = elemnumbersfromeachneighbour[n][pos+j];
                    mymaptothisdomain[n][i][recnum] = subsininterfaces[n][i][j];
                }
            }
            pos += num;
        }
    }
}

void dtracker::mapoverlapinterfaces(void)
{
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int numneighbours = myneighbours.size();
    
    int numinterfacetype = els->countinterfaceelementtypes();
    
    // Cells are guaranteed to be listed in the same order in the inner overlap region of this domain and in the outer
    // overlap region of each neighbour even if they have been renumbered during the disjoint region range definition.
    std::vector<std::vector<int>> sublistforneighbours(numneighbours), sublistfromneighbours(numneighbours);
    std::vector<std::vector<std::vector<int>>> sublistouteroverlapinterface(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        // Get the sublist on the inner overlap interface:
        std::vector<std::vector<int>> sublists(numinterfacetype, std::vector<int>(0));
        for (int i = 0; i < numinterfacetype; i++)
        {
            std::vector<std::vector<std::vector<int>>*> inneroverlapinterfaceelemlists(3, NULL);
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = myinneroverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                    inneroverlapinterfaceelemlists[dim] = prs->get(cr)->getelementlist();
            }
            std::vector<std::vector<int>>* inneroverlapelems = prs->get(myinneroverlaps[cn])->getelementlist();
        
            els->follow(inneroverlapelems, i, sublists[i], inneroverlapinterfaceelemlists);
        }
        sublistforneighbours[n] = gentools::concatenate(sublists);
        sublistforneighbours[n].resize(sublistforneighbours[n].size() + numinterfacetype);
        for (int i = 0; i < numinterfacetype; i++)
            sublistforneighbours[n][sublistforneighbours[n].size()-numinterfacetype+i] = els->count(i);
        
        // Preallocate the receive buffer based on the outer overlap interface:
        int totlen = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            std::vector<std::vector<std::vector<int>>*> outeroverlapinterfaceelems(3, NULL);
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = myouteroverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                    outeroverlapinterfaceelems[dim] = prs->get(cr)->getelementlist();
            }
            std::vector<std::vector<int>>* outeroverlapelems = prs->get(myouteroverlaps[cn])->getelementlist();
        
            els->follow(outeroverlapelems, i, sublistouteroverlapinterface[n][i], outeroverlapinterfaceelems);
            
            totlen += sublistouteroverlapinterface[n][i].size();
        }
        sublistfromneighbours[n].resize(totlen + numinterfacetype);
    }
    
    slmpi::exchange(myneighbours, sublistforneighbours, sublistfromneighbours);
    
    mymaptothisdomain = std::vector<std::vector<std::vector<int>>>(numneighbours, std::vector<std::vector<int>>(numinterfacetype, std::vector<int>(0)));
    
    // Receive the inner overlap interface of each neighbour on the outer overlap interface of this domain:
    for (int n = 0; n < numneighbours; n++)
    {
        // Preallocate the map:
        for (int i = 0; i < numinterfacetype; i++)
            mymaptothisdomain[n][i] = std::vector<int>(sublistfromneighbours[n][sublistfromneighbours[n].size()-numinterfacetype+i], -1);
            
        int pos = 0;
        for (int i = 0; i < numinterfacetype; i++)
        {
            int numsubs = sublistouteroverlapinterface[n][i].size();
            for (int j = 0; j < numsubs; j++)
                mymaptothisdomain[n][i][sublistfromneighbours[n][pos+j]] = sublistouteroverlapinterface[n][i][j];
        
            pos += numsubs;
        }
    }
}

void dtracker::createglobalnodenumbersnooverlap(void)
{
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    int numneighbours = myneighbours.size();
    
    int numnodes = nds->count();
    
    // Curvature nodes have a -1 global numbering:
    std::vector<bool> isitacornernode = els->iscornernode();
    int numcornernodes = gentools::counttrue(isitacornernode);
    
    // Number of own nodes (nodes shared with lower ranks are not own nodes):
    int numownnodes = numcornernodes;
    
    // Owner rank for each node:
    std::vector<int> nodeownerrank(numnodes, rank);
    std::vector<std::vector<int>> nodesinnoi(numneighbours, std::vector<int>(0));
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<std::vector<std::vector<int>>*> noielems(3, NULL);
        for (int d = 0; d < 3; d++)
        {
            int cr = mynooverlapinterfaces[3*cn+d];
            if (cr >= 0)
                noielems[d] = prs->get(cr)->getelementlist();
        }
        
        std::vector<bool> isinelemlists;
        int numinellist = els->istypeinelementlists(0, noielems, isinelemlists, false);
        gentools::find(isinelemlists, numinellist, nodesinnoi[n]);
        
        for (int i = 0; i < numinellist; i++)
        {
            int curnode = nodesinnoi[n][i];
            if (cn < nodeownerrank[curnode])
            {
                nodeownerrank[curnode] = cn;
                numownnodes--;
            }
        }
    }
    
    // Share with all ranks the number of own nodes:
    std::vector<int> numownnodesinthisrank = {numownnodes};
    std::vector<int> mynumownnodesinall;
    slmpi::allgather(numownnodesinthisrank, mynumownnodesinall);
    
    // Get the global number of the first own node:
    std::vector<long long int> offsetforeachrank(numranks, 0);
    for (int i = 1; i < numranks; i++)
        offsetforeachrank[i] = offsetforeachrank[i-1] + mynumownnodesinall[i-1];
    
    // Preallocate the global node numbers with a -1 default value:
    myglobalnodenumbers = std::vector<long long int>(numnodes, -1);
    
    // Make the global node numbers correct for the own nodes:

    long long int ownnodenumber = offsetforeachrank[rank];
    for (int i = 0; i < numnodes; i++)
    {
        if (isitacornernode[i] && nodeownerrank[i] == rank)
        {
            myglobalnodenumbers[i] = ownnodenumber;
            ownnodenumber++;
        }
    }
    
    // Make the global node numbers correct for the remaining nodes:
    
    std::vector<std::vector<int>> globalnodenumbersforeachneighbour(numneighbours, std::vector<int>(0));
    std::vector<std::vector<int>> globalnodenumbersfromeachneighbour(numneighbours, std::vector<int>(0));
    
    // Transfer the own nodes to the higher rank neighbours:
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
    
        // Count the number of global node numbers to send to the current neighbour:
        int numtosend = 0;
        for (int i = 0; i < nodesinnoi[n].size(); i++)
        {
            if (nodeownerrank[nodesinnoi[n][i]] == rank)
                numtosend++;
        }
        
        globalnodenumbersforeachneighbour[n].resize(2*numtosend);
        int ni = 0;
        for (int i = 0; i < nodesinnoi[n].size(); i++)
        {
            if (nodeownerrank[nodesinnoi[n][i]] == rank)
            {
                int curnode = nodesinnoi[n][i];
                globalnodenumbersforeachneighbour[n][2*ni+0] = curnode;
                globalnodenumbersforeachneighbour[n][2*ni+1] = myglobalnodenumbers[curnode]-offsetforeachrank[rank]; // referenced to offset
                ni++;
            }
        }
        
        // Count the number of global nodes to receive:
        int numtoreceive = 0;
        for (int i = 0; i < nodesinnoi[n].size(); i++)
        {
            if (nodeownerrank[nodesinnoi[n][i]] == cn)
                numtoreceive++;
        }
        globalnodenumbersfromeachneighbour[n].resize(2*numtoreceive);
    }
    
    slmpi::exchange(myneighbours, globalnodenumbersforeachneighbour, globalnodenumbersfromeachneighbour);
    
    // Set the missing global node numbers:
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < globalnodenumbersfromeachneighbour[n].size()/2; i++)
            myglobalnodenumbers[mymaptothisdomain[n][0][globalnodenumbersfromeachneighbour[n][2*i+0]]] = globalnodenumbersfromeachneighbour[n][2*i+1]+offsetforeachrank[myneighbours[n]];
    }
}

void dtracker::createglobalnodenumbersoverlap(void)
{
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    int numneighbours = myneighbours.size();
    
    int numnodes = nds->count();
    
    // Split neighbours into lower ranks and higher ranks (neighbours are sorted ascendingly):
    std::vector<int> lowerranks, higherranks;
    gentools::split(myneighbours, rank, lowerranks, higherranks);
    
    // Curvature nodes have a -1 global numbering:
    std::vector<bool> isitacornernode = els->iscornernode();
    int numcornernodes = gentools::counttrue(isitacornernode);
    
    // Number of own nodes (nodes shared with lower ranks are not own nodes):
    int numownnodes = numcornernodes;
    
    // Owner rank for each node:
    std::vector<int> nodeownerrank(numnodes, rank);
    
    // Nodes in the no-overlap interfaces, inner overlaps and outer overlaps:
    std::vector<std::vector<int>> nodesinnoiseenfromio(numneighbours, std::vector<int>(0));
    std::vector<std::vector<int>> nodesinnoiseenfromoo(numneighbours, std::vector<int>(0));
    std::vector<std::vector<int>> nodesinio(numneighbours, std::vector<int>(0));
    std::vector<std::vector<int>> nodesinoo(numneighbours, std::vector<int>(0));

    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        std::vector<std::vector<std::vector<int>>*> noielems(3, NULL);
        for (int d = 0; d < 3; d++)
        {
            int cr = mynooverlapinterfaces[3*cn+d];
            if (cr >= 0)
                noielems[d] = prs->get(cr)->getelementlist();
        }

        std::vector<std::vector<int>>* ioelems = prs->get(myinneroverlaps[cn])->getelementlist();
        std::vector<std::vector<int>>* ooelems = prs->get(myouteroverlaps[cn])->getelementlist();
        
        els->follow(ioelems, 0, nodesinnoiseenfromio[n], noielems);
        els->follow(ooelems, 0, nodesinnoiseenfromoo[n], noielems);
        els->follow(ioelems, 0, nodesinio[n]);
        els->follow(ooelems, 0, nodesinoo[n]);
    }
    
    // Give away the ownership of the outer overlaps:
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < nodesinoo[n].size(); i++)
        {
            int curnode = nodesinoo[n][i];
            if (nodeownerrank[curnode] == rank)
            {
                nodeownerrank[curnode] = -1;
                numownnodes--;
            }
        }
    }
    
    // Take back the ownership of the no-overlap interfaces:
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < nodesinnoiseenfromio[n].size(); i++)
        {
            int curnode = nodesinnoiseenfromio[n][i];
            if (nodeownerrank[curnode] != rank)
            {
                nodeownerrank[curnode] = rank;
                numownnodes++;
            }
        }
    }
    
    // Give the ownership of the no-overlap interface nodes to the lowest touching rank:
    for (int n = 0; n < lowerranks.size(); n++)
    {
        int cn = myneighbours[n];
        for (int i = 0; i < nodesinnoiseenfromio[n].size(); i++)
        {
            int curnode = nodesinnoiseenfromio[n][i];
            if (cn < nodeownerrank[curnode])
            {
                nodeownerrank[curnode] = cn;
                numownnodes--;
            }
        }
    }
    
    // Share with all ranks the number of own nodes:
    std::vector<int> numownnodesinthisrank = {numownnodes};
    std::vector<int> mynumownnodesinall;
    slmpi::allgather(numownnodesinthisrank, mynumownnodesinall);
    
    // Get the global number of the first own node:
    std::vector<long long int> offsetforeachrank(numranks, 0);
    for (int i = 1; i < numranks; i++)
        offsetforeachrank[i] = offsetforeachrank[i-1] + mynumownnodesinall[i-1];
    
    // Preallocate the global node numbers with a -1 default value:
    myglobalnodenumbers = std::vector<long long int>(numnodes, -1);
    
    // Make the global node numbers correct for the own nodes:

    long long int ownnodenumber = offsetforeachrank[rank];
    for (int i = 0; i < numnodes; i++)
    {
        if (isitacornernode[i] && nodeownerrank[i] == rank)
        {
            myglobalnodenumbers[i] = ownnodenumber;
            ownnodenumber++;
        }
    }

    // Make the global node numbers correct for the no-overlap interfaces:

    // Transfer the own nodes to the higher rank neighbours:
    std::vector<std::vector<int>> globalnoinodesforneighbours(numneighbours, std::vector<int>(0));
    std::vector<std::vector<int>> globalnoinodesfromneighbours(numneighbours, std::vector<int>(0));
    for (int n = 0; n < numneighbours; n++)
    {
        int numnoinodes = nodesinnoiseenfromio[n].size();
        if (n < lowerranks.size())
            globalnoinodesfromneighbours[n].resize(numnoinodes);
        else
        {
            globalnoinodesforneighbours[n].resize(numnoinodes);
            // Not owned nodes will send a -1 global number:
            for (int i = 0; i < numnoinodes; i++)
            {
                long long int gnn = myglobalnodenumbers[nodesinnoiseenfromio[n][i]];
                if (gnn >= 0)
                    globalnoinodesforneighbours[n][i] = gnn-offsetforeachrank[rank]; // referenced to offset
                else
                    globalnoinodesforneighbours[n][i] = -1;
            }
        }
    }
    
    slmpi::exchange(myneighbours, globalnoinodesforneighbours, globalnoinodesfromneighbours);
    
    for (int n = 0; n < lowerranks.size(); n++)
    {
        int cn = myneighbours[n];
        int numnoinodes = nodesinnoiseenfromoo[n].size();
        
        for (int i = 0; i < numnoinodes; i++)
        {
            int gnn = globalnoinodesfromneighbours[n][i];
            if (gnn >= 0)
                myglobalnodenumbers[nodesinnoiseenfromoo[n][i]] = gnn+offsetforeachrank[cn];
        }
    }
    
    // Make the global node numbers correct for the outer overlaps:

    // Transfer the inner overlap nodes of this domain to the outer overlap nodes of every neighbour:
    std::vector<std::vector<int>> globalnodesforneighbours(numneighbours, std::vector<int>(0));
    std::vector<std::vector<int>> globalnodesfromneighbours(numneighbours, std::vector<int>(0));
    for (int n = 0; n < numneighbours; n++)
    {
        globalnodesforneighbours[n].resize(2*nodesinio[n].size());
        for (int i = 0; i < nodesinio[n].size(); i++)
        {
            int curowner = nodeownerrank[nodesinio[n][i]];
            globalnodesforneighbours[n][2*i+0] = curowner;
            globalnodesforneighbours[n][2*i+1] = myglobalnodenumbers[nodesinio[n][i]]-offsetforeachrank[curowner]; // referenced to offset
        }
    
        globalnodesfromneighbours[n].resize(2*nodesinoo[n].size());
    }
    
    slmpi::exchange(myneighbours, globalnodesforneighbours, globalnodesfromneighbours);
    
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < nodesinoo[n].size(); i++)
            myglobalnodenumbers[nodesinoo[n][i]] = globalnodesfromneighbours[n][2*i+1]+offsetforeachrank[globalnodesfromneighbours[n][2*i+0]];
    }
}

void dtracker::setconnectivity(std::vector<int>& neighbours, std::vector<int>& nooverlapinterfaces)
{
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    int meshdim = getrawmesh()->getmeshdimension();

    if (3*neighbours.size() != nooverlapinterfaces.size())
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: expected a number of no-overlap interface regions equal to 3 x number of neighbours (one region per interface element dimension, -1 if none)" << std::endl;
        log.error();
    }

    // Make sure there is no duplicate and not the rank itself:
    int numneighbours = 0;
    myisneighbour = std::vector<bool>(numranks, false);
    mynooverlapinterfaces = std::vector<int>(3*numranks, -1);
    for (int i = 0; i < neighbours.size(); i++)
    {
        int n = neighbours[i];
        if (n >= 0 && n < numranks && myisneighbour[n] == false && n != rank)
        {
            myisneighbour[n] = true;
            
            for (int d = 0; d < 3; d++)
            {
                int ci = nooverlapinterfaces[3*i+d];
                
                if (ci >= 0)
                {
                    prs->errorundefined({ci});
                    int elemdim = prs->get(ci)->getelementdimension();
                    
                    if (elemdim >= meshdim)
                    {
                        logs log;
                        log.msg() << "Error in 'dtracker' object: provided an interface physical region with " << elemdim << "D elements (this is not an interface)" << std::endl;
                        log.error();
                    }
                    if (elemdim != d)
                    {
                        logs log;
                        log.msg() << "Error in 'dtracker' object: expected a physical region with " << d << "D elements at index 3*" << i << "+" << d << std::endl;
                        log.error();
                    }
                    
                    mynooverlapinterfaces[3*n+d] = ci;
                }
            }
            
            numneighbours++;
        }
        else
        {
            logs log;
            log.msg() << "Error in 'dtracker' object: neighbours provided must be unique numbers between 0 and " << numranks << " and cannot include the domain rank itself" << std::endl;
            log.error();
        }
    }
    
    myneighbours = std::vector<int>(numneighbours);
    
    int index = 0;
    for (int i = 0; i < numranks; i++)
    {
        if (myisneighbour[i])
        {
            myneighbours[index] = i;
            index++;
        }
    }
}

void dtracker::discoverconnectivity(int numtrialelements, int verbosity)
{        
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int rank = slmpi::getrank();
    int numranks = slmpi::count();

    // Define the interface of the domain with its neighbours:
    int skinregion = prs->getmaxphysicalregionnumber()+1;
    int wholeneighbourinterface = skinregion+1;

    regiondefiner regdef(*nds, *els, *prs);
    regdef.regionskin(skinregion, -1);
    // This domain might not be in contact with the global geometry skin:
    if (myglobalgeometryskin >= 0)
        regdef.regionexclusion(wholeneighbourinterface, skinregion, {myglobalgeometryskin});
    else
        wholeneighbourinterface = skinregion;
    regdef.defineregions();
    
    std::vector<std::vector<int>> interfaceelems = *(prs->get(wholeneighbourinterface)->getelementlist());    
    // Remove the construction regions:
    prs->remove({skinregion, wholeneighbourinterface}, false);

    // Avoid any conflict between DDM-created regions and the other regions:
    std::vector<int> lastpr = {prs->getmaxphysicalregionnumber()};
    slmpi::max(lastpr);
    int firstnewpr = lastpr[0]+1;
    
    int meshdim = getrawmesh()->getmeshdimension();
    numtrialelements = std::max(1,numtrialelements);
    
    std::vector<bool> isnodeininterface, isedgeininterface;
    int numnodesininterface = els->istypeinelementlists(0, {&interfaceelems}, isnodeininterface, false);
    int numedgesininterface = els->istypeinelementlists(1, {&interfaceelems}, isedgeininterface, false);
    std::vector<int> interfacenodelist, interfaceedgelist;
    gentools::find(isnodeininterface, numnodesininterface, interfacenodelist);
    gentools::find(isedgeininterface, numedgesininterface, interfaceedgelist);

    std::vector<physicalregion*> physregsvec(3*numranks, NULL);

    std::vector<int> allnei = {};

    int numits = 0;
    while (true)
    {   
        std::vector<double> elembarys;
        els->getbarycenters(&interfaceelems, elembarys);
        
        std::vector<int> neighboursfound;
        std::vector<int> allnumelementsininterface = discoversomeneighbours(numtrialelements, elembarys, neighboursfound);

        if (allnumelementsininterface.size() == 0)
            break;
            
        if (allnumelementsininterface == allnei)
        {
            logs log;
            log.msg() << "Error in 'dtracker' object: connectivity discovery algorithm failed because some interface elements could not be found on any other domain" << std::endl;
            log.error();
        }
        allnei = allnumelementsininterface;
        
        // Discover the cell-1 dimension interfaces:
        std::vector<int> inneighbour;
        discoverinterfaces(neighboursfound, elembarys, allnumelementsininterface, inneighbour);
        
        // Add the matched elements to their physical region and remove them from the list:
        int index = 0;
        for (int i = 0; i < 8; i++)
        {
            int numelemstokeep = 0;
            for (int j = 0; j < interfaceelems[i].size(); j++)
            {
                int elem = interfaceelems[i][j];
                int curneighbour = inneighbour[index];
                
                if (curneighbour != -1)
                {
                    if (physregsvec[3*curneighbour+(meshdim-1)] == NULL)
                    {
                        physregsvec[3*curneighbour+(meshdim-1)] = prs->get(firstnewpr);
                        firstnewpr++;   
                    }
                        
                    physregsvec[3*curneighbour+(meshdim-1)]->addelement(i, elem);
                }
                else
                {
                    interfaceelems[i][numelemstokeep] = elem;
                    numelemstokeep++;
                }
                index++; 
            }
            interfaceelems[i].resize(numelemstokeep);
        }

        numits++;
    }
    
    
    // Discover the neighbours touching with lower dimension interfaces (none in 1D, only nodes in 2D, nodes and edges in 3D):
    
    // One entry per rank (any rank is a neighbour candidate):
    std::vector<std::vector<bool>> isnodeinneighbours(numranks, std::vector<bool>(0));
    std::vector<std::vector<bool>> isedgeinneighbours(numranks, std::vector<bool>(0));
    for (int r = 0; r < numranks; r++)
    {
        if (physregsvec[3*r+(meshdim-1)] != NULL)
        {
            els->istypeinelementlists(0, {physregsvec[3*r+(meshdim-1)]->getelementlist()}, isnodeinneighbours[r], false);
            els->istypeinelementlists(1, {physregsvec[3*r+(meshdim-1)]->getelementlist()}, isedgeinneighbours[r], false);
        }
    }
    
    std::vector<std::vector<bool>> wasnodeinneighbours = isnodeinneighbours;
    std::vector<std::vector<bool>> wasedgeinneighbours = isedgeinneighbours;
        
    // Discover all cross-interfaces:
    int numcrossits = 0;
    while (true)
    {
        bool isanyfound = discovercrossinterfaces(interfacenodelist, interfaceedgelist, isnodeinneighbours, isedgeinneighbours);
        if (isanyfound == false)
            break;
        numcrossits++;
    }

    // Add first the edges then the nodes that are not in the edges:
    for (int r = 0; r < numranks; r++)
    {
        for (int i = 0; i < isedgeinneighbours[r].size(); i++)
        {
            if (isedgeinneighbours[r][i] == true && (wasedgeinneighbours[r].size() == 0 || wasedgeinneighbours[r][i] == false))
            {
                int nodea = els->getsubelement(0, 1, i, 0);
                int nodeb = els->getsubelement(0, 1, i, 1);
                
                // The nodes should not be added to the interface region since they are part of an added edge:
                isnodeinneighbours[r][nodea] = false;
                isnodeinneighbours[r][nodeb] = false;
                
                if (physregsvec[3*r+1] == NULL)
                {
                    physregsvec[3*r+1] = prs->get(firstnewpr);
                    firstnewpr++;   
                }
                    
                physregsvec[3*r+1]->addelement(1, i);
            }
        }

        for (int i = 0; i < isnodeinneighbours[r].size(); i++)
        {
            if (isnodeinneighbours[r][i] == true && (wasnodeinneighbours[r].size() == 0 || wasnodeinneighbours[r][i] == false))
            {
                if (physregsvec[3*r+0] == NULL)
                {
                    physregsvec[3*r+0] = prs->get(firstnewpr);
                    firstnewpr++;
                }
                    
                physregsvec[3*r+0]->addelement(0, i);
            }
        }
    }
    
    // Create connectivity containers:
    myneighbours = {};
    myisneighbour = std::vector<bool>(numranks, false);
    mynooverlapinterfaces = std::vector<int>(3*numranks, -1);
    for (int r = 0; r < numranks; r++)
    {
        if (physregsvec[3*r+0] != NULL || physregsvec[3*r+1] != NULL || physregsvec[3*r+2] != NULL)
        {
            if (r == rank)
            {
                logs log;
                log.msg() << "Error in 'dtracker' object: rank " << r << " is a neighbour of itself (this is unexpected)" << std::endl;
                log.error();
            }
        
            myneighbours.push_back(r);
            myisneighbour[r] = true;
            for (int dim = 0; dim < 3; dim++)
            {
                if (physregsvec[3*r+dim] != NULL)
                    mynooverlapinterfaces[3*r+dim] = physregsvec[3*r+dim]->getnumber();
            }
        }
    }
    
    if (verbosity > 0)
        std::cout << "Found all neighbours with " << numits << " set" << gentools::getplurals(numits) << " of " << numtrialelements << " trial element" << gentools::getplurals(numtrialelements) << " and " << numcrossits << " propagation step" << gentools::getplurals(numcrossits) << std::endl;
}

void dtracker::overlap(void)
{
    if (isoverlap())
    {   
        defineinneroverlaps();
        exchangeoverlaps();
        exchangephysicalregions();
        defineouteroverlapinterfaces();
        defineinneroverlapinterfaces();
    }
}

void dtracker::mapinterfaces(void)
{
    if (isoverlap())
        mapoverlapinterfaces();
    else
        mapnooverlapinterfaces();
}

void dtracker::createglobalnodenumbers(void)
{
    if (isoverlap())
        createglobalnodenumbersoverlap();
    else
        createglobalnodenumbersnooverlap();
}


int dtracker::countneighbours(void)
{
    return myneighbours.size();
}

std::vector<int> dtracker::getneighbours(void)
{
    return myneighbours;
}

int dtracker::getneighbour(int neighbourindex)
{
    if (neighbourindex >= 0 && neighbourindex < myneighbours.size())
        return myneighbours[neighbourindex];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: asked for a neighbour at an index larger than the number of neighbours" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

bool dtracker::isneighbour(int neighbour)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size())
        return myisneighbour[neighbour];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: asked if rank " << neighbour << " is a neighbour of rank " << slmpi::getrank() << " but there are only " << slmpi::count() << " ranks in total" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getnooverlapinterface(int neighbour, int elementdimension)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size() && elementdimension >= 0 && elementdimension < 3)
        return mynooverlapinterfaces[3*neighbour+elementdimension];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: requested on rank " << slmpi::getrank() << " the " << elementdimension << "D no-overlap interface to neighbour rank " << neighbour << " but there are only " << slmpi::count() << " ranks in total" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getinneroverlapinterface(int neighbour, int elementdimension)
{
    if (neighbour >= 0 && 3*neighbour < myinneroverlapinterfaces.size() && elementdimension >= 0 && elementdimension < 3)
        return myinneroverlapinterfaces[3*neighbour+elementdimension];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested inner overlap interface region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getouteroverlapinterface(int neighbour, int elementdimension)
{
    if (neighbour >= 0 && 3*neighbour < myouteroverlapinterfaces.size() && elementdimension >= 0 && elementdimension < 3)
        return myouteroverlapinterfaces[3*neighbour+elementdimension];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested outer overlap interface region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getinneroverlap(int neighbour)
{
    if (neighbour >= 0 && neighbour < myinneroverlaps.size())
        return myinneroverlaps[neighbour];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested inner overlap region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getouteroverlap(int neighbour)
{
    if (neighbour >= 0 && neighbour < myouteroverlaps.size())
        return myouteroverlaps[neighbour];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested outer overlap region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getinneroverlapskin(int neighbour)
{
    if (neighbour >= 0 && neighbour < myinneroverlapskins.size())
        return myinneroverlapskins[neighbour];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested inner overlap skin region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

int dtracker::getouteroverlapskin(int neighbour)
{
    if (neighbour >= 0 && neighbour < myouteroverlapskins.size())
        return myouteroverlapskins[neighbour];
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested outer overlap skin region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> dtracker::getnooverlapinterface(int neighbour)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size())
    {
        std::vector<int> output = {};
        for (int d = 0; d < 3; d++)
        {
            int cpr = mynooverlapinterfaces[3*neighbour+d];
            if (cpr >= 0)
                output.push_back(cpr);
        }
        return output;
    }
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: requested on rank " << slmpi::getrank() << " the no-overlap interface to neighbour rank " << neighbour << " but there are only " << slmpi::count() << " ranks in total" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> dtracker::getinneroverlapinterface(int neighbour)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size())
    {
        std::vector<int> output = {};
        for (int d = 0; d < 3; d++)
        {
            int cpr = myinneroverlapinterfaces[3*neighbour+d];
            if (cpr >= 0)
                output.push_back(cpr);
        }
        return output;
    }
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested inner overlap interface region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> dtracker::getouteroverlapinterface(int neighbour)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size())
    {
        std::vector<int> output = {};
        for (int d = 0; d < 3; d++)
        {
            int cpr = myouteroverlapinterfaces[3*neighbour+d];
            if (cpr >= 0)
                output.push_back(cpr);
        }
        return output;
    }
    else
    {
        logs log;
        log.msg() << "Error in 'dtracker' object: cannot provide the requested outer overlap interface region" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::vector<std::vector<int>>>* dtracker::getmap(void)
{
    return &mymaptothisdomain;
}

long long int* dtracker::getglobalnodenumbers(void)
{
    return myglobalnodenumbers.data();
}

void dtracker::writeglobalnodenumbers(std::string filename)
{
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();

    int numnodes = nds->count();

    // Convert to int:
    std::vector<int> nodenums(numnodes);
    for (int i = 0; i < numnodes; i++)
        nodenums[i] = myglobalnodenumbers[i];

    els->write(filename, 0, gentools::getequallyspaced(0, 1, numnodes), nodenums);
}

std::vector<int> dtracker::listddmregions(void)
{
    std::vector<int> prlist = {};
    
    int numneighbours = myneighbours.size();
    
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = myneighbours[n];
        
        if (mynooverlapinterfaces.size() > 0)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = mynooverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                    prlist.push_back(cr);
            }
        }
        
        if (myinneroverlaps.size() > 0 && myinneroverlaps[cn] >= 0)
            prlist.push_back(myinneroverlaps[cn]);
        if (myouteroverlaps.size() > 0 && myouteroverlaps[cn] >= 0)
            prlist.push_back(myouteroverlaps[cn]);
        if (myinneroverlapskins.size() > 0 && myinneroverlapskins[cn] >= 0)
            prlist.push_back(myinneroverlapskins[cn]);
        if (myouteroverlapskins.size() > 0 && myouteroverlapskins[cn] >= 0)
            prlist.push_back(myouteroverlapskins[cn]);
            
        if (myinneroverlapinterfaces.size() > 0)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = myinneroverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                    prlist.push_back(cr);
            }
        }
        if (myouteroverlapinterfaces.size() > 0)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = myouteroverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                    prlist.push_back(cr);
            }
        }
    }
    
    return prlist;
}

std::vector<bool> dtracker::isddmdisjointregion(void)
{
    disjointregions* drs = getrawmesh()->getdisjointregions();
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    std::vector<bool> output(drs->count(), false);

    std::vector<int> prlist = listddmregions();
    for (int i = 0; i < prlist.size(); i++)
    {
        std::vector<bool> isindr = prs->get(prlist[i])->getdefinition();
        for (int d = 0; d < output.size(); d++)
            output[d] = (output[d] || isindr[d]);
    }
    return output;
}

std::vector<bool> dtracker::isdisjointregioninnooverlap(void)
{
    disjointregions* drs = getrawmesh()->getdisjointregions();
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    std::vector<bool> isinnooverlap(drs->count(), true);

    if (isoverlap())
    {
        for (int n = 0; n < myneighbours.size(); n++)
        {
            std::vector<bool> drsdef = prs->get(myouteroverlaps[myneighbours[n]])->getdefinition();
            for (int d = 0; d < drsdef.size(); d++)
            {
                if (drsdef[d])
                    isinnooverlap[d] = false;
            }
        }
        for (int n = 0; n < myneighbours.size(); n++)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = mynooverlapinterfaces[3*myneighbours[n]+dim];
                if (cr >= 0)
                {
                    std::vector<bool> drsdef = prs->get(cr)->getdefinition();
                    for (int d = 0; d < drsdef.size(); d++)
                    {
                        if (drsdef[d])
                            isinnooverlap[d] = true;
                    }
                }
            }
        }
    }

    return isinnooverlap;
}

std::vector<bool> dtracker::isdisjointregionowned(void)
{
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    
    std::vector<bool> isowned = isdisjointregioninnooverlap();
    
    for (int n = 0; n < myneighbours.size(); n++)
    {
        int cn = myneighbours[n];
        if (cn < rank)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                int cr = mynooverlapinterfaces[3*cn+dim];
                if (cr >= 0)
                {
                    std::vector<bool> drsdef = prs->get(cr)->getdefinition();
                    for (int d = 0; d < drsdef.size(); d++)
                    {
                        if (drsdef[d])
                            isowned[d] = false;
                    }
                }
            }
        }
    }

    return isowned;
}

void dtracker::print(void)
{
    std::vector<int> lens = {(int)myneighbours.size()};
    std::vector<int> alllens;
    slmpi::gather(0, lens, alllens);
    std::vector<int> allneighbours;
    slmpi::gather(0, myneighbours, allneighbours, alllens);

    if (slmpi::getrank() == 0)
    {
        int index = 0;
        for (int r = 0; r < myisneighbour.size(); r++)
        {
            std::cout << "Rank " << r << " has " << alllens[r] << " neighbours";
            for (int n = 0; n < alllens[r]; n++)
                std::cout << " " << allneighbours[index+n];
            std::cout << std::endl;
                
            index += alllens[r];
        }
    }
}

void dtracker::writeinterfaces(std::string filename)
{
    if (filename.size() >= 5)
    {
        // Get the extension:
        std::string fileext = filename.substr(filename.size()-4,4);
        // Get the file name without the extension:
        std::string noext = filename.substr(0, filename.size()-4);

        int rank = slmpi::getrank();
        
        for (int r = 0; r < myisneighbour.size(); r++)
        {
            if (myisneighbour[r])
            {
                std::string cfn = noext + "dom" + std::to_string(rank) + "nei" + std::to_string(r);
                for (int dim = 0; dim < 3; dim++)
                {
                    if (mynooverlapinterfaces[3*r+dim] != -1)
                        expression(rank).write(mynooverlapinterfaces[3*r+dim], cfn + "dim" + std::to_string(dim) + fileext, 1);
                }
            }
        }

        return;
    }
    
    logs log;
    log.msg() << "Error in 'dtracker' object: cannot write to file '" << filename << "'." << std::endl;
    log.msg() << "Supported output formats are .vtk (ParaView), .vtu (ParaView) and .pos (GMSH)." << std::endl;
    log.error();
}

