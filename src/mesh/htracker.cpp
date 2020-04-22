#include "htracker.h"
#include "myalgorithm.h"
#include "lagrangeformfunction.h"


htracker::htracker(std::vector<int> numelemspertype)
{
    originalcount = numelemspertype;
    
    for (int i = 0; i < 8; i++)
        numleaves += originalcount[i];
    
    // Initialize the tree:
    splitdata = std::vector<bool>(numleaves, false);
}

int htracker::countleaves(void)
{
    return numleaves;
}

int htracker::getmaxdepth(void)
{
    return maxdepth;
}

void htracker::resetcursor(void)
{
    cursorposition = 0;
    currentdepth = 0;
    curtypeorigcountindex = 0;
    // Find the first not empty type:
    parenttypes = std::vector<int>(maxdepth+1,0);
    while (originalcount[parenttypes[0]] == 0)
        parenttypes[0]++;
    indexesinclusters = std::vector<int>(maxdepth+1,0);
}

int htracker::next(void)
{
    int throughedgenum = -1;

    // If not split:
    if (splitdata[cursorposition] == false)
    {
        cursorposition++;
        
        if (currentdepth > 0)
        {
            indexesinclusters[currentdepth]++;
            // For pyramids the subelement type might change from 4 to 7:
            if (parenttypes[currentdepth-1] == 7 && indexesinclusters[currentdepth] == 4)
                parenttypes[currentdepth] = 7;
                
            while (currentdepth > 0 && indexesinclusters[currentdepth] == numsubelems[parenttypes[currentdepth-1]])
                currentdepth--;
        }
        if (currentdepth == 0)
        {
            curtypeorigcountindex++;
            // Skip empty types:
            while (curtypeorigcountindex == originalcount[parenttypes[0]])
            {
                parenttypes[0]++;
                curtypeorigcountindex = 0;
            }
        }
    }
    else
    {
        // For pyramids the subelement type might change from 4 to 7:
        if (currentdepth > 0 && parenttypes[currentdepth-1] == 7 && indexesinclusters[currentdepth] == 4)
            parenttypes[currentdepth] = 7;    
        // Split pyramids give 4 tets and 6 pyramids:
        if (parenttypes[currentdepth] != 7)
            parenttypes[currentdepth+1] = parenttypes[currentdepth];
        else
            parenttypes[currentdepth+1] = 4;

        if (currentdepth > 0)
            indexesinclusters[currentdepth]++;
        indexesinclusters[currentdepth+1] = 0;

        cursorposition++;
        if (parenttypes[currentdepth] == 4)
        {
            throughedgenum = 0;
            if (splitdata[cursorposition])
                throughedgenum += 2;
            if (splitdata[cursorposition+1])
                throughedgenum += 1;
                
            cursorposition += 2;
        }
        
        currentdepth++;
    }
    
    return throughedgenum;
}

bool htracker::isatleaf(void)
{
    return (splitdata[cursorposition] == 0);
}

int htracker::countsplits(void)
{
    return currentdepth;
}

int htracker::gettype(void)
{
    return parenttypes[currentdepth];
}

int htracker::getparenttype(void)
{
    if (currentdepth > 0)
        return parenttypes[currentdepth-1];
    else
        return -1;
}

int htracker::getindexincluster(void)
{
    return indexesinclusters[currentdepth];
}

void htracker::countsplits(std::vector<int>& numsplits)
{
    numsplits = std::vector<int>(numleaves);
    
    resetcursor();
    
    for (int i = 0; i < numleaves; i++)
    {    
        // Move to next leaf:
        while (not(isatleaf()))
            next();
            
        numsplits[i] = currentdepth;
        
        next();
    }
}

void htracker::gettype(std::vector<int>& types)
{
    types = std::vector<int>(numleaves);
    
    resetcursor();
    
    for (int i = 0; i < numleaves; i++)
    {    
        // Move to next leaf:
        while (not(isatleaf()))
            next();
            
        types[i] = parenttypes[currentdepth];
        
        next();
    }
}

std::vector<int> htracker::countintypes(void)
{
    std::vector<int> output(8,0);
    
    resetcursor();
    
    for (int i = 0; i < numleaves; i++)
    {    
        // Move to next leaf:
        while (not(isatleaf()))
            next();
            
        output[parenttypes[currentdepth]]++;
        
        next();
    }
    
    return output;
}

void htracker::fix(std::vector<int>& operations)
{
    // Number of grouping requests in a cluster:
    std::vector<int> ngr(maxdepth+1);
    
    resetcursor();
    
    int ln = -1; // leaf number
    
    while (true)
    {
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
        
        // Reset cluster info:
        if (ic == 0)
            ngr[ns] = 0;
        
        if (isatleaf())
        {
            ln++;
            
            if (operations[ln] == -1)
            {
                operations[ln] = 0;
                ngr[ns]++;
            }
            
            // Group the cluster if all have requested it:
            if (ns > 0 && ngr[ns] == numsubelems[parenttypes[ns-1]])
            {
                // In this situation all leaves are consecutive:
                for (int i = 0; i < numsubelems[parenttypes[ns-1]]; i++)
                    operations[ln-i] = -1;
            }
        }
        
        if (ln == numleaves-1)
            break;
        
        next();
    }
}

void htracker::adapt(std::vector<int>& operations, std::vector<int>& throughedgenums)
{
    // Calculate an upper bound for the size of the new 'splitdata' vector:
    int upperbound = splitdata.size();
    for (int i = 0; i < operations.size(); i++)
    {
        if (operations[i] == 1)
            upperbound += 10; // 10 is the max size increase in a split
    }
    std::vector<bool> newsplitdata(upperbound, false);

    // Number of grouping requests in a cluster:
    std::vector<int> ngr(maxdepth+1);
    
    resetcursor();
    
    int newnumleaves = numleaves;
    int newmaxdepth = 0;
    
    int ln = -1; // leaf number
    int ni = 0; // newsplitdata index
    int cte = -1; // current throughedge
    
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
        
        // Reset cluster info:
        if (ic == 0)
            ngr[ns] = 0;

        // Write the throughedge number:
        if (cte != -1)
        {
            if (cte%2 == 1)
                newsplitdata[ni+1] = true;
            if (((cte-cte%2)/2)%2 == 1)
                newsplitdata[ni] = true;
            ni += 2;
        }
        
        if (not(isatleaf()))
        {
            newsplitdata[ni] = true;
            ni++;
        }
        else
        {
            ln++;
            
            // Split:
            if (operations[ln] == 1)
            {
                newnumleaves += numsubelems[t]-1;
                
                newsplitdata[ni] = true;
                ni++;
                if (t == 4)
                {
                    if (throughedgenums[ln]%2 == 1)
                        newsplitdata[ni+1] = true;
                    if (((throughedgenums[ln]-throughedgenums[ln]%2)/2)%2 == 1)
                        newsplitdata[ni] = true;
                    ni += 2;
                }
                // Nothing to change for subelements (initialized at all false + group puts it to false anyways):
                ni += numsubelems[t];
                
                if (newmaxdepth < ns+1)
                    newmaxdepth = ns+1;
            }
            
            if (operations[ln] == -1)
                ngr[ns]++;
                
            // Group/unchanged:
            if (operations[ln] < 1)
            {
                newsplitdata[ni] = false;
                ni++;
            }
            
            // Group the cluster if all have requested it:
            if (ns > 0 && ic == numsubelems[parenttypes[ns-1]]-1)
            {
                int md = ns;
                if (ngr[ns] == numsubelems[parenttypes[ns-1]])
                {
                    md--;
                    int pt = parenttypes[ns-1];
                    newnumleaves -= numsubelems[pt]-1;
                    ni -= numsubelems[pt];
                    if (pt == 4)
                        ni -= 2;
                    newsplitdata[ni-1] = false;
                }
                if (newmaxdepth < md)
                    newmaxdepth = md;
            }
        }
        
        if (ln == numleaves-1)
            break;
        
        cte = next();
    }
    
    newsplitdata.resize(ni);
    splitdata = newsplitdata;
    numleaves = newnumleaves;
    maxdepth = newmaxdepth;
}

void htracker::print(void)
{
    std::cout << "#leaves is " << numleaves << " | max #splits is " << maxdepth << std::endl;

    for (int i = 0; i < splitdata.size(); i++)
        std::cout << splitdata[i];
        
    std::cout << std::endl;
}

int htracker::countbits(void)
{
    return splitdata.size();
}

void htracker::countsons(std::vector<int>& numsons)
{
    int num = 0;
    for (int i = 0; i < 8; i++)
        num += originalcount[i];
        
    numsons = std::vector<int>(8*num,0);
    
    resetcursor();
    
    int ln = -1;
    int index = -1;
    while (true)
    { 
        if (currentdepth == 0)
            index++;
            
        if (isatleaf())
        {
            ln++;
            numsons[8*index+parenttypes[currentdepth]]++;
        }
        
        if (ln == numleaves-1)
            break;
        
        next();
    }
}

void htracker::getadapted(int curvatureorder, std::vector<std::vector<double>>& oc, std::vector<std::vector<double>>& arc, std::vector<std::vector<double>>& ac)
{
    std::vector<int> nit = countintypes();
    std::vector<int> nn(8); // number of corner nodes
    std::vector<int> ncn(8); // number of corner nodes

    // Preallocate:
    std::vector<element> els(8);
    arc = std::vector<std::vector<double>>(8, std::vector<double>(0));
    ac = std::vector<std::vector<double>>(8, std::vector<double>(0));
    std::vector<std::vector<double>> crc(8); // corner ref coords
    for (int i = 0; i < 8; i++)
    {
        els[i] = element(i);
        nn[i] = els[i].countnodes();
        element curvedel(i,curvatureorder);
        ncn[i] = curvedel.countcurvednodes();
        
        arc[i] = std::vector<double>(3*nn[i]*nit[i]);
        ac[i] = std::vector<double>(3*nn[i]*nit[i]);
        lagrangeformfunction lff(i,1,{});
        crc[i] = lff.getnodecoordinates();
    }

    // parentrc[depth][indexincluster]:
    std::vector<std::vector<std::vector<double>>> parentrc(maxdepth+1, std::vector<std::vector<double>>(10));
    
    std::vector<double> originalelemcoords;
    
    resetcursor();
    
    int ln = -1;
    std::vector<int> iarc(8,0); // indexes in arc
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        if (ns == 0)
        {
            parentrc[0] = {crc[t]};
            originalelemcoords = std::vector<double>(3*nn[t]);
            for (int i = 0; i < 3*nn[t]; i++)
                originalelemcoords[i] = oc[t][curtypeorigcountindex*3*ncn[t]+i];
        }
    
        if (isatleaf())
        {
            ln++;
            
            std::vector<double> realcoords;
            if (ns == 0)
                realcoords = originalelemcoords;
            else
                realcoords = els[parenttypes[0]].calculatecoordinates(parentrc[ns][ic], originalelemcoords);
            
            for (int i = 0; i < parentrc[ns][ic].size(); i++)
            {
                arc[t][iarc[t]+i] = parentrc[ns][ic][i];
                ac[t][iarc[t]+i] = realcoords[i];
            }
            iarc[t] += parentrc[ns][ic].size();
            
        
            if (ln == numleaves-1)
                break;
            
            next();
        }
        else
        {
            int throughedgenum = next();
                    
            std::vector<std::vector<double>> cornerrefcoords;
            els[t].fullsplit(cornerrefcoords, throughedgenum);
            
            int ind = 0;
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < cornerrefcoords[i].size()/nn[i]/3; j++)
                {
                    std::vector<double> cc(3*nn[i]);
                    for (int k = 0; k < cc.size(); k++)
                        cc[k] = cornerrefcoords[i][j*3*nn[i]+k];
                    
                    parentrc[ns+1][ind] = els[t].calculatecoordinates(cc, parentrc[ns][ic]);

                    ind++;
                }
            }
        }
    }
}

void htracker::getadaptedcoordinates(int curvatureorder, std::vector<std::vector<double>>& oc, std::vector<std::vector<double>>& ac, std::vector<double> noisethreshold)
{
    std::vector<int> nn(8);
    std::vector<int> ncn(8);
    std::vector<int> ne(8);
    std::vector<element> straightelements(8);
    std::vector<element> curvedelements(8);
    std::vector<std::vector<double>> curvedcoords(8);
    for (int i = 0; i < 8; i++)
    {
        straightelements[i] = element(i);
        curvedelements[i] = element(i,curvatureorder);
        nn[i] = curvedelements[i].countnodes();
        ncn[i] = curvedelements[i].countcurvednodes();
        ne[i] = curvedelements[i].countedges();
        lagrangeformfunction lff(i,curvatureorder,{});
        curvedcoords[i] = lff.getnodecoordinates();
    }

    // Get the reference ('arc') and real ('ac') element corner coordinates after all fullsplit adaptation:
    std::vector<std::vector<double>> cornerac;
    std::vector<std::vector<double>> cornerarc;
    getadapted(curvatureorder, oc, cornerarc, cornerac);
    
    // Compute the barycenter coordinates of all nodes and edges (first edges then nodes):
    int numnodes = 0, numedges = 0;
    for (int i = 0; i < 8; i++)
    {
        numnodes += cornerac[i].size()/3;
        numedges += ne[i]*cornerac[i].size()/nn[i]/3;
    }
    std::vector<double> barys(3*numedges+3*numnodes);

    int ce = 0, cn = 0;
    for (int i = 0; i < 8; i++)
    {
        element el(i);
        std::vector<int> edgenodedef = el.getedgesdefinitionsbasedonnodes();
    
        int num = cornerac[i].size()/nn[i]/3;
        for (int j = 0; j < num; j++)
        {
            for (int e = 0; e < ne[i]; e++)
            {
                int na = edgenodedef[2*e+0];
                int nb = edgenodedef[2*e+1];
                
                barys[3*ce+0] = 0.5*(cornerac[i][3*nn[i]*j+3*na+0] + cornerac[i][3*nn[i]*j+3*nb+0]);
                barys[3*ce+1] = 0.5*(cornerac[i][3*nn[i]*j+3*na+1] + cornerac[i][3*nn[i]*j+3*nb+1]);
                barys[3*ce+2] = 0.5*(cornerac[i][3*nn[i]*j+3*na+2] + cornerac[i][3*nn[i]*j+3*nb+2]);
                
                ce++;
            }
            for (int e = 0; e < nn[i]; e++)
            {
                barys[3*numedges+3*cn+0] = cornerac[i][3*nn[i]*j+3*e+0];
                barys[3*numedges+3*cn+1] = cornerac[i][3*nn[i]*j+3*e+1];
                barys[3*numedges+3*cn+2] = cornerac[i][3*nn[i]*j+3*e+2];
                
                cn++;
            }
        }
    }
    
    
    // Sort the barycenter coordinates:
    std::vector<double> sortedbarys(barys.size());
    std::vector<int> reorderingvector;
    myalgorithm::stablecoordinatesort(noisethreshold, barys, reorderingvector);
    for (int i = 0; i < reorderingvector.size(); i++)
    {
        sortedbarys[3*i+0] = barys[3*reorderingvector[i]+0];
        sortedbarys[3*i+1] = barys[3*reorderingvector[i]+1];
        sortedbarys[3*i+2] = barys[3*reorderingvector[i]+2];
    }
    std::vector<int> renum(reorderingvector.size());
    for (int i = 0; i < reorderingvector.size(); i++)
        renum[reorderingvector[i]] = i;
    
    // Remove duplicated barycenters:
    std::vector<int> renumberingvector;
    int numunique = myalgorithm::removeduplicatedcoordinates(noisethreshold, sortedbarys, renumberingvector);
    
    
    // Assign a unique edge number for each edge:
    std::vector<int> edgenumbers(numedges);
    for (int i = 0; i < numedges; i++)
        edgenumbers[i] = renumberingvector[renum[i]];
    
    // Calculate which edges must be split:
    std::vector<bool> isanodeatnum(numedges+numnodes,false);
    for (int i = 0; i < numnodes; i++)
        isanodeatnum[renumberingvector[renum[numedges+i]]] = true;

    std::vector<bool> isedgesplit(numedges,false);
    for (int i = 0; i < numedges; i++)
    {
        if (isanodeatnum[edgenumbers[i]])
            isedgesplit[i] = true;
    }
    
    
    // Split the transition elements:
    std::vector<int> numsons;
    countsons(numsons);
    
    ac = std::vector<std::vector<double>>(8, std::vector<double>(0));
    // No-transition size:
    std::vector<int> nts(8,0);
    for (int i = 0; i < 8; i++)
        nts[i] = ncn[i] * cornerarc[i].size()/nn[i];
    // Preallocate to an upper bound:
    ac[1] = std::vector<double>(nts[1]);
    ac[2] = std::vector<double>(4*nts[2]+3*nts[3]); // NO: ASK FOR ELEMENT OBJECT TO GIVE THAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ac[3] = std::vector<double>(4*nts[3]);
    ac[4] = std::vector<double>(8*nts[4]);
    // Undefined for other elements for now:
    ac[5] = std::vector<double>(0);
    ac[6] = std::vector<double>(0);
    ac[7] = std::vector<double>(0);
        
    
    int origelemindex = 0;
    std::vector<int> indexincoords(8,0); // index in working element
    std::vector<int> indexintransitioncoords(8,0); // index in working element
    std::vector<int> firstedge(8,0); // first edge in working element
    for (int i = 0; i < 7; i++)
        firstedge[i+1] = firstedge[i] + ne[i] * cornerarc[i].size()/nn[i]/3;
    // Loop on all original elements:
    for (int i = 0; i < 8; i++)
    {
        for (int oe = 0; oe < originalcount[i]; oe++)
        {
            // Get the coordinates of the current original element:
            std::vector<double> origcoords(3*ncn[i]);
            for (int j = 0; j < 3*ncn[i]; j++)
                origcoords[j] = oc[i][3*ncn[i]*oe+j];
        
            // Loop on all elements types in the fullsplit-subelements of the original element:
            for (int j = 0; j < 8; j++)
            {
                for (int e = 0; e < numsons[8*origelemindex+j]; e++)
                {
                    // Get the corner ref. coords. of the current subelement:
                    std::vector<double> currefcoords(3*nn[j]);
                    for (int k = 0; k < 3*nn[j]; k++)
                        currefcoords[k] = cornerarc[j][indexincoords[j]+k];
                    indexincoords[j] += 3*nn[j];
                    
                    // Get the edge numbers and edge splits for the current subelement:
                    std::vector<int> curedgenums(ne[j]);
                    std::vector<bool> curisedgesplit(ne[j]);
                    for (int k = 0; k < ne[j]; k++)
                    {
                        curedgenums[k] = edgenumbers[firstedge[j]+k];
                        curisedgesplit[k] = isedgesplit[firstedge[j]+k];
                    }
                    
                    int splitnum = myalgorithm::binarytoint(curisedgesplit);
                    std::vector<std::vector<int>> splitrefnums = straightelements[j].split(splitnum, curedgenums);
                  
                    // Loop on all subelements in the transition element:
                    for (int si = 0; si < 8; si++)
                    {
                        std::vector<double> splitrefcoords;
                        straightelements[j].numstorefcoords(splitrefnums[si], splitrefcoords);
                    
                        for (int se = 0; se < splitrefcoords.size()/nn[si]/3; se++)
                        {
                            // Get the ref. coords. of the current transition-subelement:
                            std::vector<double> curcoords(3*nn[si]);
                            for (int k = 0; k < 3*nn[si]; k++)
                                curcoords[k] = splitrefcoords[se*nn[si]*3+k];
                    
                            // Bring inside the untransitioned element (if split at all):
                            if (splitnum == 0)
                                curcoords = currefcoords;
                            else
                                curcoords = straightelements[j].calculatecoordinates(curcoords, currefcoords);
                            
                            // Make curved:
                            if (curvatureorder > 1)
                                curcoords = straightelements[si].calculatecoordinates(curvedcoords[si], curcoords);
                                
                            // Calculate actual coordinates: 
                            curcoords = curvedelements[i].calculatecoordinates(curcoords, origcoords);

                            for (int k = 0; k < curcoords.size(); k++)
                                ac[si][indexintransitioncoords[si]+k] = curcoords[k];
                                
                            indexintransitioncoords[si] += curcoords.size();
                        }
                    }
                    firstedge[j] += ne[j];
                }
            }
            origelemindex++;
        }
    }
    
    // Fit to size:
    for (int i = 0; i < 8; i++)
        ac[i].resize(indexintransitioncoords[i]);
    
}

