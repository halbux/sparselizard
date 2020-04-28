#include "htracker.h"
#include "myalgorithm.h"
#include "lagrangeformfunction.h"


htracker::htracker(int curvatureorder, std::vector<int> numelemspertype)
{
    originalcurvatureorder = curvatureorder;
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

void htracker::getoriginalelementnumber(std::vector<int>& oen)
{
    oen = std::vector<int>(numleaves);
    
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
            oen[ln] = index;
        }
        
        if (ln == numleaves-1)
            break;
        
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

void htracker::adapt(std::vector<int>& operations)
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
                    newsplitdata[ni] = true;
                    newsplitdata[ni+1] = true;
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

void htracker::getadapted(std::vector<std::vector<double>>& oc, std::vector<std::vector<double>>& arc, std::vector<std::vector<double>>& ac)
{
    std::vector<int> nit = countintypes();
    std::vector<int> nn(8); // number of corner nodes
    std::vector<int> ncn(8); // number of curved nodes

    // Preallocate:
    std::vector<element> els(8);
    std::vector<element> curvedels(8);
    arc = std::vector<std::vector<double>>(8, std::vector<double>(0));
    std::vector<std::vector<double>> crc(8); // corner ref coords
    std::vector<std::vector<double>> curvedrc(8);
    for (int i = 0; i < 8; i++)
    {
        els[i] = element(i);
        nn[i] = els[i].countnodes();
        curvedels[i] = element(i,originalcurvatureorder);
        ncn[i] = curvedels[i].countcurvednodes();
        
        arc[i] = std::vector<double>(3*nn[i]*nit[i]);
        lagrangeformfunction lff(i,1,{});
        crc[i] = lff.getnodecoordinates();
        lagrangeformfunction lffc(i,originalcurvatureorder,{});
        curvedrc[i] = lffc.getnodecoordinates();
    }
    ac = arc;

    // Corner reference coordinates - parentrc[depth][indexincluster]:
    std::vector<std::vector<std::vector<double>>> parentrc(maxdepth+1, std::vector<std::vector<double>>(10));
    
    resetcursor();
    
    int ln = -1; // leaf number
    std::vector<int> iarc(8,0); // indexes in arc
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ot = parenttypes[0];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        if (ns == 0)
            parentrc[0] = {crc[t]};
        
        if (isatleaf())
        {
            ln++;
            
            std::vector<double> physcoords = els[parenttypes[0]].calculatecoordinates(parentrc[ns][ic], oc[parenttypes[0]], curtypeorigcountindex*3*ncn[parenttypes[0]], ns == 0);
            
            for (int i = 0; i < parentrc[ns][ic].size(); i++)
            {
                arc[t][iarc[t]+i] = parentrc[ns][ic][i];
                ac[t][iarc[t]+i] = physcoords[i];
            }
            
            iarc[t] += parentrc[ns][ic].size();
        
            if (ln == numleaves-1)
                break;
            
            next();
        }
        else
        {
            int throughedgenum = next();
            
            // Calculate best through-edge if not yet defined:
            if (t == 4 && throughedgenum == 3)
            {
                std::vector<double> coordsforten = els[t].calculatecoordinates(curvedrc[t], parentrc[ns][ic], originalcurvatureorder == 1);
                coordsforten = curvedels[ot].calculatecoordinates(coordsforten, oc[ot], curtypeorigcountindex*ncn[ot]*3, ns == 0);
                
                throughedgenum = curvedels[t].choosethroughedge(coordsforten);
                
                // Add it to 'splitdata' once and for all:
                if (throughedgenum%2 == 1)
                    splitdata[cursorposition-1] = true;
                if (((throughedgenum-throughedgenum%2)/2)%2 == 1)
                    splitdata[cursorposition-2] = true;
            }
                    
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

void htracker::getadaptedcoordinates(std::vector<std::vector<double>>& oc, std::vector<std::vector<double>>& ac, std::vector<std::vector<int>>& leafnums, std::vector<double> noisethreshold)
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
        curvedelements[i] = element(i,originalcurvatureorder);
        nn[i] = curvedelements[i].countnodes();
        ncn[i] = curvedelements[i].countcurvednodes();
        ne[i] = curvedelements[i].countedges();
        lagrangeformfunction lff(i,originalcurvatureorder,{});
        curvedcoords[i] = lff.getnodecoordinates();
    }


    // Get the reference ('arc') and physical ('ac') element corner coordinates after all fullsplit adaptation:
    std::vector<std::vector<double>> cornerac;
    std::vector<std::vector<double>> cornerarc;
    getadapted(oc, cornerarc, cornerac);
    
    
    // Assign unique edge numbers and deduce edge splits:
    std::vector<int> edgenumbers;
    std:vector<bool> isedgesplit;
    myalgorithm::assignedgenumbers(cornerac, edgenumbers, isedgesplit, noisethreshold);
    

    // Preallocate output containers:
    ac = std::vector<std::vector<double>>(8, std::vector<double>(0));
    teorc = std::vector<std::vector<double>>(8, std::vector<double>(0));
    transitionelemsleafnums = std::vector<std::vector<int>>(8, std::vector<int>(0));
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
    
    for (int i = 0; i < 8; i++)
    {
        transitionelemsleafnums[i] = std::vector<int>(ac[i].size()/ncn[i]/3);
        teorc[i] = std::vector<double>(nn[i]*ac[i].size()/ncn[i]);
    }
    
    resetcursor();
    
    int ln = -1; // leaf number
    std::vector<int> iarc(8,0); // indexes in arc
    std::vector<int> firstedge(8,0); // first edge in working element
    for (int i = 0; i < 7; i++)
        firstedge[i+1] = firstedge[i] + ne[i] * cornerarc[i].size()/nn[i]/3;
    std::vector<int> iac(8,0); // index in ac
    while (true)
    {
        while (not(isatleaf()))
            next();
    
        ln++;
     
        int t = parenttypes[currentdepth];
        int ot = parenttypes[0];
        
        // Get the edge numbers and edge splits for the current subelement:
        std::vector<int> curedgenums(ne[t]);
        std::vector<bool> curisedgesplit(ne[t]);
        for (int i = 0; i < ne[t]; i++)
        {
            curedgenums[i] = edgenumbers[firstedge[t]+i];
            curisedgesplit[i] = isedgesplit[firstedge[t]+i];
        }
        
        int splitnum = myalgorithm::binarytoint(curisedgesplit);
        std::vector<std::vector<int>> splitrefnums = straightelements[t].split(splitnum, curedgenums);
      
        // Loop on all transition elements:
        for (int si = 0; si < 8; si++)
        {
            if (splitrefnums[si].size() == 0)
                continue;
        
            std::vector<double> splitrefcoords;
            straightelements[t].numstorefcoords(splitrefnums[si], splitrefcoords);
        
            for (int se = 0; se < splitrefcoords.size()/nn[si]/3; se++)
            {
                // Get the ref. coords. of the current transition element:
                std::vector<double> curcoords(3*nn[si]);
                for (int i = 0; i < 3*nn[si]; i++)
                    curcoords[i] = splitrefcoords[se*nn[si]*3+i];
        
                // Bring inside the untransitioned element (if split at all):
                curcoords = straightelements[t].calculatecoordinates(curcoords, cornerarc[t], iarc[t], splitnum == 0);
                
                for (int i = 0; i < curcoords.size(); i++)
                    teorc[si][nn[si]*iac[si]/ncn[si]+i] = curcoords[i];
                
                // Make curved:
                if (originalcurvatureorder > 1)
                    curcoords = straightelements[si].calculatecoordinates(curvedcoords[si], curcoords);
                    
                // Calculate actual coordinates: 
                curcoords = curvedelements[ot].calculatecoordinates(curcoords, oc[ot], 3*ncn[ot]*curtypeorigcountindex);

                for (int i = 0; i < curcoords.size(); i++)
                    ac[si][iac[si]+i] = curcoords[i];
                    
                transitionelemsleafnums[si][iac[si]/ncn[si]/3] = ln;
                
                iac[si] += 3*ncn[si];
                
            }
        }
        
        firstedge[t] += ne[t];
        iarc[t] += 3*nn[t];
    
        if (ln == numleaves-1)
            break;
        
        next();
    }
    
    // Fit to size:
    for (int i = 0; i < 8; i++)
    {
        if (iac[i] > 0)
        {
            ac[i].resize(iac[i]);
            teorc[i].resize(nn[i]*iac[i]/ncn[i]);
            transitionelemsleafnums[i].resize(iac[i]/ncn[i]/3);
        }
        // ELSE RESIZE TO 0 OR ALREADY EMPTY?
    }
    leafnums = transitionelemsleafnums;
}

void htracker::inoriginal(std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& oad, std::vector<double>& orc)
{
    std::vector<int> oen;
    getoriginalelementnumber(oen);

    // Preallocate:
    int numrc = 0, no = 0;
    std::vector<int> nn(8);
    std::vector<element> straightels(8);
    for (int i = 0; i < 8; i++)
    {
        straightels[i] = element(i);
        nn[i] = straightels[i].countnodes();
        numrc += rc[i].size()/3;
        no += originalcount[i]; 
    }   
    oad = std::vector<int>(no+1,0);
    orc = std::vector<double>(3*numrc);
    
    // Create the 'oad' address vector:
    std::vector<int> cnt(no,0);
    for (int i = 0; i < 8; i++)
    {
        // This is needed because .size() gives an unsigned int --> .size()-1 underflows
        if (ad[i].size() == 0)
            continue;
        for (int j = 0; j < ad[i].size()-1; j++)
        {
            int origelem = oen[transitionelemsleafnums[i][j]];
            cnt[origelem] += ad[i][j+1]-ad[i][j];
        }
    }
    
    for (int i = 1; i <= no; i++)
        oad[i] = oad[i-1]+cnt[i-1];
        
    // Loop on each transition element:
    std::vector<int> index(no,0);
    for (int i = 0; i < 8; i++)
    {
        // This is needed because .size() gives an unsigned int --> .size()-1 underflows
        if (ad[i].size() == 0)
            continue;
        for (int j = 0; j < ad[i].size()-1; j++)
        {
            int nr = (ad[i][j+1]-ad[i][j])/3;
            
            std::vector<double> currefs(3*nr);       
            for (int k = 0; k < 3*nr; k++)
                currefs[k] = rc[i][ad[i][j]+k];
                
            currefs = straightels[i].calculatecoordinates(currefs, teorc[i], 3*nn[i]*j);
    
            int origelem = oen[transitionelemsleafnums[i][j]];
            for (int k = 0; k < 3*nr; k++)
                orc[oad[origelem]+index[origelem]+k] = currefs[k];
       
            index[origelem] += 3*nr;
        }
    }
}

void htracker::tostorage(void)
{
    
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

