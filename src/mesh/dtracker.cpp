#include "dtracker.h"
#include "nodes.h"
#include "elements.h"
#include "physicalregions.h"


dtracker::dtracker(std::shared_ptr<rawmesh> rm)
{
    if (rm == NULL)
    {
        std::cout << "Error in 'dtracker' object: cannot provide a NULL rawmesh pointer" << std::endl;
        abort();
    }
    
    myrawmesh = rm;
}

std::shared_ptr<rawmesh> dtracker::getrawmesh(void)
{
    if (myrawmesh.expired() == false)
        return myrawmesh.lock();
    else
    {
        std::cout << "Error in 'dtracker' object: cannot get rawmesh object (weak pointer is expired)" << std::endl;
        abort();
    }
}

std::vector<int> dtracker::discoversomeneighbours(int numtrialelements, std::vector<double>& interfaceelembarys, std::vector<int>& neighboursfound)
{ 
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    neighboursfound = {};
    
    int numelementsininterface = interfaceelembarys.size()/3;
    
    std::vector<double> trialbarys(3*numtrialelements, 0.0); // must have this size in any case, barycenters for not alive ranks will be removed below
    if (numelementsininterface > 0)
        myalgorithm::pickcandidates(numtrialelements, interfaceelembarys, trialbarys); // ok if not unique
    
    // Push this in the mpi call:
    trialbarys.push_back(myalgorithm::exactinttodouble(std::min(numelementsininterface,1)));
    
    std::vector<double> alltrialbarys;
    slmpi::allgather(trialbarys, alltrialbarys);
    
    std::vector<double> allisalive = myalgorithm::extract(alltrialbarys, 3*numtrialelements+1, 3*numtrialelements);
    
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
        myalgorithm::findcoordinates(interfaceelembarys, alltrialbarys, isfound);
        
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
    
    std::vector<int> allnumelementsininterface = myalgorithm::extract(allneighbourlist, numtrialelements+1, numtrialelements);
    
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
    int numelementsininterface = interfaceelembarys.size()/3;
    
    inneighbour = std::vector<int>(numelementsininterface, -1);
    
    int numneighbours = neighbours.size();
    if (numneighbours == 0)
        return;
    
    // Place in one vector all received candidate coordinates:
    int totnumcand = 0;
    for (int n = 0; n < numneighbours; n++)
        totnumcand += allnumelementsininterface[neighbours[n]];
        
    std::vector<double> candidatebarys(3*totnumcand);
    
    // Exchange all interface barycenter coordinates with every neighbour:
    std::vector<int> receivelens(numneighbours);
    std::vector<double*> receivebuffers(numneighbours);
    
    int pos = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int len = 3*allnumelementsininterface[neighbours[n]];
        receivelens[n] = len;
        receivebuffers[n] = &candidatebarys[pos];
        pos += len;
    }
    slmpi::exchange(neighbours, 3*numelementsininterface, &interfaceelembarys[0], receivelens, receivebuffers);
    
    // Find matches:
    std::vector<int> posfound;
    myalgorithm::findcoordinates(interfaceelembarys, candidatebarys, posfound);
    
    pos = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int len = allnumelementsininterface[neighbours[n]];
        for (int i = 0; i < len; i++)
        {
            int pf = posfound[pos+i];
            if (pf != -1)
                inneighbour[pf] = neighbours[n];
        }
        pos += len;
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
            neighbours.push_back(i);
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
    
    // Package edge and node barycenters to both neighbours in each neighbour pair:
    std::vector<std::vector<double>> packaged(numneighbours*numneighbours, std::vector<double>(0));
    for (int i = 0; i < numneighbours*numneighbours; i++)
    {
        int totnumtosend = numnodesinneighbourpair[i] + numedgesinneighbourpair[i];
        
        packaged[i] = std::vector<double>(3*totnumtosend+2);
        packaged[i][0] = 3*myalgorithm::exactinttodouble(numnodesinneighbourpair[i]);
        packaged[i][1] = 3*myalgorithm::exactinttodouble(numedgesinneighbourpair[i]);
        
        myalgorithm::selectcoordinates(nodesinneighbourpair[i], *ncs, &(packaged[i][2]));
        myalgorithm::selectcoordinates(edgesinneighbourpair[i], *edgebarys, &(packaged[i][2+3*numnodesinneighbourpair[i]]));
    }
    
    std::vector<std::vector<double>> dataforeachneighbour;
    myalgorithm::pack(neighbours, packaged, dataforeachneighbour);
    
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
        sendbuffers[n] = &(dataforeachneighbour[n][0]);
    }
    slmpi::exchange(neighbours, sendlens, receivelens);
    
    for (int n = 0; n < numneighbours; n++)
    {
        datafromeachneighbour[n].resize(receivelens[n]);
        receivebuffers[n] = &(datafromeachneighbour[n][0]);
    }
    slmpi::exchange(neighbours, sendlens, sendbuffers, receivelens, receivebuffers);
    
    // Unpack:
    std::vector< std::vector<std::vector<double>> > unpackedcoords(numneighbours);
    std::vector<std::vector<int>> candidateneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
        candidateneighbours[n] = myalgorithm::unpack(datafromeachneighbour[n], unpackedcoords[n]);
    
    std::vector<double> allreceivednodecoords, allreceivededgecoords;
    std::vector<std::vector<int>> numnodesingroup, numedgesingroup;
    myalgorithm::split(unpackedcoords, allreceivednodecoords, allreceivededgecoords, numnodesingroup, numedgesingroup);
    
    std::vector<double> interfacenodesbarys;
    els->getbarycenters(0, interfacenodelist, interfacenodesbarys);
    std::vector<double> interfaceedgesbarys;
    els->getbarycenters(1, interfaceedgelist, interfaceedgesbarys);
    
    std::vector<int> posnodefound;
    int nnf = myalgorithm::findcoordinates(interfacenodesbarys, allreceivednodecoords, posnodefound);
    std::vector<int> posedgefound;
    int nef = myalgorithm::findcoordinates(interfaceedgesbarys, allreceivededgecoords, posedgefound);

    // All nodes and edges must have been found or something went wrong:
    if (nnf != posnodefound.size() || nef != posedgefound.size())
    {
        std::cout << "Error in 'dtracker' object: algorithm to detect cross-interfaces did not succeed (at least one node or edge barycenter was not matched on a target neighbour)" << std::endl;
        abort();
    }
    
    // Update 'isnodeinneighbours' and 'isedgeinneighbours':
    int nodeindex = 0;
    int edgeindex = 0;
    
    int isanynewadded = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        for (int i = 0; i < unpackedcoords[n].size(); i++)
        {
            if (unpackedcoords[n][i].size() == 0)
                continue;
                
            int nn = numnodesingroup[n][i]/3;
            int ne = numedgesingroup[n][i]/3;
            
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

void dtracker::setconnectivity(std::vector<int>& neighbours, std::vector<int>& nooverlapinterfaces)
{
    physicalregions* prs = getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    int meshdim = getrawmesh()->getmeshdimension();

    if (3*neighbours.size() != nooverlapinterfaces.size())
    {
        std::cout << "Error in 'dtracker' object: expected a number of no-overlap interface regions equal to 3 x number of neighbours (one region per interface element dimension, -1 if none)" << std::endl;
        abort();
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
                        std::cout << "Error in 'dtracker' object: provided an interface physical region with " << elemdim << "D elements (this is not an interface)" << std::endl;
                        abort();
                    }
                    if (elemdim != d)
                    {
                        std::cout << "Error in 'dtracker' object: expected a physical region with " << d << "D elements at index 3*" << i << "+" << d << std::endl;
                        abort();
                    }
                    
                    mynooverlapinterfaces[3*n+d] = ci;
                }
            }
            
            numneighbours++;
        }
        else
        {
            std::cout << "Error in 'dtracker' object: neighbours provided must be unique numbers between 0 and " << numranks << " and cannot include the domain rank itself" << std::endl;
            abort();
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

void dtracker::discoverconnectivity(int nooverlapinterface, int numtrialelements, int verbosity)
{        
    nodes* nds = getrawmesh()->getnodes();
    elements* els = getrawmesh()->getelements();
    physicalregions* prs = getrawmesh()->getphysicalregions();

    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    int meshdim = getrawmesh()->getmeshdimension();
    
    std::vector<std::vector<int>> interfaceelems = *(prs->get(nooverlapinterface)->getelementlist());

    std::vector<bool> isnodeininterface, isedgeininterface;
    int numnodesininterface = els->istypeinelementlist(0, &interfaceelems, isnodeininterface);
    int numedgesininterface = els->istypeinelementlist(1, &interfaceelems, isedgeininterface);
    std::vector<int> interfacenodelist, interfaceedgelist;
    myalgorithm::find(isnodeininterface, numnodesininterface, interfacenodelist);
    myalgorithm::find(isedgeininterface, numedgesininterface, interfaceedgelist);

    std::vector<physicalregion*> physregsvec(3*numranks, NULL);

    std::vector<int> allnei = {};

    int numits = 0;
    while (true)
    {   
        std::vector<double> elembarys;
        els->getbarycenters(&interfaceelems, elembarys);
        
        std::vector<int> neighboursfound;
        std::vector<int> allnumelementsininterface = discoversomeneighbours(std::max(1,numtrialelements), elembarys, neighboursfound);

        int numneighboursfound = neighboursfound.size();
        
        if (allnumelementsininterface.size() == 0)
            break;
            
        if (allnumelementsininterface == allnei)
        {
            std::cout << "Error in 'dtracker' object: connectivity discovery algorithm failed because some interface elements could not be found on any other domain" << std::endl;
            abort();
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
                        physregsvec[3*curneighbour+(meshdim-1)] = prs->get(prs->getmaxphysicalregionnumber()+1);
                        
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
            els->istypeinelementlist(0, physregsvec[3*r+(meshdim-1)]->getelementlist(), isnodeinneighbours[r]);
            els->istypeinelementlist(1, physregsvec[3*r+(meshdim-1)]->getelementlist(), isedgeinneighbours[r]);
        }
    }
    
    int numcrossits = 0;
    while (true)
    {
        std::vector<std::vector<bool>> wasnodeinneighbours = isnodeinneighbours;
        std::vector<std::vector<bool>> wasedgeinneighbours = isedgeinneighbours;

        bool isanyfound = discovercrossinterfaces(interfacenodelist, interfaceedgelist, isnodeinneighbours, isedgeinneighbours);

        if (isanyfound == false)
            break;

        // First add the newly discovered edges:
        int numnodes = nds->count();
        int numedges = els->count(1);
        for (int r = 0; r < numranks; r++)
        {
            if (r == rank || isnodeinneighbours.size() == 0 && isedgeinneighbours.size() == 0)
                continue;
            
            // Nodes and edges might be added:
            if (wasnodeinneighbours[r].size() == 0)
                wasnodeinneighbours[r] = std::vector<bool>(numnodes, false);
            if (wasedgeinneighbours[r].size() == 0)
                wasedgeinneighbours[r] = std::vector<bool>(numedges, false);
                
            for (int i = 0; i < isedgeinneighbours[r].size(); i++)
            {
                // A new edge must be added to the neighbour interface:
                if (isedgeinneighbours[r][i] == true && wasedgeinneighbours[r][i] == false)
                {
                    int nodea = els->getsubelement(0, 1, i, 0);
                    int nodeb = els->getsubelement(0, 1, i, 1);
                    
                    isnodeinneighbours[r][nodea] = true;
                    isnodeinneighbours[r][nodeb] = true;
                    // The node should not be added to the interface region since it is part of an added edge:
                    wasnodeinneighbours[r][nodea] = true;
                    wasnodeinneighbours[r][nodeb] = true;
                    
                    // Add the edge to the physical region:
                    if (physregsvec[3*r+1] == NULL)
                        physregsvec[3*r+1] = prs->get(prs->getmaxphysicalregionnumber()+1);
                        
                    physregsvec[3*r+1]->addelement(1, i);
                }
            }

            for (int i = 0; i < isnodeinneighbours[r].size(); i++)
            {
                // A new node must be added to the neighbour interface:
                if (isnodeinneighbours[r][i] == true && wasnodeinneighbours[r][i] == false)
                {
                    // Add the node to the physical region:
                    if (physregsvec[3*r+0] == NULL)
                        physregsvec[3*r+0] = prs->get(prs->getmaxphysicalregionnumber()+1);
                        
                    physregsvec[3*r+0]->addelement(0, i);
                }
            }
        }
        
        numcrossits++;
    }
    
    // Create connectivity containers:
    myneighbours = {};
    myisneighbour = std::vector<bool>(numranks, false);
    mynooverlapinterfaces = std::vector<int>(3*numranks, -1);
    for (int r = 0; r < numranks; r++)
    {
        if (physregsvec[3*r+0] != NULL || physregsvec[3*r+1] != NULL || physregsvec[3*r+2] != NULL)
        {
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
        std::cout << "Found all neighbour domains with " << numits << " set" << myalgorithm::getplurals(numits) << " of " << numtrialelements << " trial element" << myalgorithm::getplurals(numtrialelements) << " and " << numcrossits << " propagation step" << myalgorithm::getplurals(numcrossits) << std::endl;
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
        std::cout << "Error in 'dtracker' object: asked for a neighbour at an index larger than the number of neighbours" << std::endl;
        abort();
    }
}

bool dtracker::isneighbour(int neighbour)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size())
        return myisneighbour[neighbour];
    else
    {
        std::cout << "Error in 'dtracker' object: asked if rank " << neighbour << " is a neighbour of rank " << slmpi::getrank() << " but there are only " << slmpi::count() << " ranks in total" << std::endl;
        abort();
    }
}

int dtracker::getnooverlapinterface(int neighbour, int elementdimension)
{
    if (neighbour >= 0 && neighbour < myisneighbour.size() && elementdimension >= 0 && elementdimension < 3)
        return mynooverlapinterfaces[3*neighbour+elementdimension];
    else
    {
        std::cout << "Error in 'dtracker' object: requested on rank " << slmpi::getrank() << " the " << elementdimension << "D no-overlap interface to neighbour rank " << neighbour << " but there are only " << slmpi::count() << " ranks in total" << std::endl;
        abort();
    }
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
    
    std::cout << "Error in 'dtracker' object: cannot write to file '" << filename << "'." << std::endl;
    std::cout << "Supported output formats are .vtk (ParaView), .vtu (ParaView) and .pos (GMSH)." << std::endl;
    abort();
}

