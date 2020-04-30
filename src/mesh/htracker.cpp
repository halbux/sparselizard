#include "htracker.h"
#include "myalgorithm.h"
#include "lagrangeformfunction.h"


htracker::htracker(elements* origelems, int curvatureorder, std::vector<int> numelemspertype)
{
    myoriginalelements = origelems;
    
    originalcurvatureorder = curvatureorder;
    originalcount = numelemspertype;
    
    if (curvatureorder == -1)
        originalcurvatureorder = myoriginalelements->getcurvatureorder();
    if (numelemspertype.size() == 0)
    {
        int dim = myoriginalelements->getdimension();
        originalcount = std::vector<int>(8,0);
        for (int i = 0; i < 8; i++)
        {
            element myelem(i);
            if (myelem.getelementdimension() == dim)
                originalcount[i] = myoriginalelements->count(i);
        }
    }
        
    myelems = std::vector<element>(8);
    mycurvedelems = std::vector<element>(8);
    nn = std::vector<int>(8);
    ncn = std::vector<int>(8);
    straightrefcoords = std::vector<std::vector<double>>(8);
    for (int i = 0; i < 8; i++)
    {
        numleaves += originalcount[i];
        lagrangeformfunction lff(i,1,{});
        straightrefcoords[i] = lff.getnodecoordinates();   
        myelems[i] = element(i);
        mycurvedelems[i] = element(i,originalcurvatureorder);
        nn[i] = myelems[i].countnodes();
        ncn[i] = mycurvedelems[i].countcurvednodes();
    }
    
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

void htracker::resetcursor(bool calcrefcoords)
{
    isrefcalc = calcrefcoords;

    cursorposition = 0;
    currentdepth = 0;
    origindexintype = 0;
    // Find the first not empty type:
    parenttypes = std::vector<int>(maxdepth+1,0);
    while (originalcount[parenttypes[0]] == 0)
        parenttypes[0]++;
    indexesinclusters = std::vector<int>(maxdepth+1,0);
    parentrefcoords = std::vector<std::vector<std::vector<double>>>(maxdepth+1, std::vector<std::vector<double>>(10));
}

int htracker::next(void)
{
    int throughedgenum = -1;
    
    if (isrefcalc && currentdepth == 0)
        parentrefcoords[0] = {straightrefcoords[parenttypes[0]]};

    // If not split we are at a leaf:
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
            origindexintype++;
            // Skip empty types:
            while (origindexintype == originalcount[parenttypes[0]])
            {
                parenttypes[0]++;
                origindexintype = 0;
            }
        }
    }
    else
    {
        cursorposition++;
        // Get the through-edge number (if any):
        if (parenttypes[currentdepth] == 4)
        {
            throughedgenum = 0;
            if (splitdata[cursorposition])
                throughedgenum += 2;
            if (splitdata[cursorposition+1])
                throughedgenum += 1;
                
            // If it is 3 we have to define it:
            if (throughedgenum == 3)
            {
                std::vector<double> origelemcoords = myoriginalelements->getnodecoordinates(parenttypes[0], origindexintype);
                throughedgenum = mycurvedelems[4].choosethroughedge(origelemcoords);
                // Write it to the tree:
                if (throughedgenum%2 == 1)
                    splitdata[cursorposition+1] = true;
                if (((throughedgenum-throughedgenum%2)/2)%2 == 1)
                    splitdata[cursorposition] = true;
            }
            cursorposition += 2;
        }
        
        // Optionally compute the current reference coordinates in the original element:
        if (isrefcalc)
        {
            std::vector<std::vector<double>> cornerrefcoords;
            myelems[parenttypes[currentdepth]].fullsplit(cornerrefcoords, throughedgenum);
            
            int ind = 0;
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < cornerrefcoords[i].size()/nn[i]/3; j++)
                {
                    std::vector<double> cc(3*nn[i]);
                    for (int k = 0; k < cc.size(); k++)
                        cc[k] = cornerrefcoords[i][j*3*nn[i]+k];
                    
                    parentrefcoords[currentdepth+1][ind] = myelems[parenttypes[currentdepth]].calculatecoordinates(cc, parentrefcoords[currentdepth][indexesinclusters[currentdepth]]);

                    ind++;
                }
            }
        }
    
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

std::vector<double> htracker::getreferencecoordinates(void)
{
    if (currentdepth == 0)
        return straightrefcoords[parenttypes[0]];
    else
        return parentrefcoords[currentdepth][indexesinclusters[currentdepth]];
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

        // Mark the throughedge number as undefined (value 3):
        if (cte != -1)
        {
            newsplitdata[ni+1] = true;
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

void htracker::atleaves(std::vector<std::vector<double>>& arc, std::vector<std::vector<double>>& apc, bool withphysicals)
{
    std::vector<int> nit = countintypes();

    // Preallocate:
    arc = std::vector<std::vector<double>>(8, std::vector<double>(0));
    for (int i = 0; i < 8; i++)
        arc[i] = std::vector<double>(3*nn[i]*nit[i]);
    if (withphysicals)
        apc = arc;

    resetcursor(true);
    
    std::vector<double> oc;
    
    int ln = -1; // leaf number
    std::vector<int> iarc(8,0); // indexes in arc
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ot = parenttypes[0];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        if (ns == 0 && ic == 0)
            oc = myoriginalelements->getnodecoordinates(t, origindexintype);
    
        if (isatleaf())
        {
            ln++;
            
            std::vector<double> refcoords = getreferencecoordinates();
            std::vector<double> physcoords;
            if (withphysicals)
                physcoords = myelems[parenttypes[0]].calculatecoordinates(refcoords, oc, 0, ns == 0);
            
            for (int i = 0; i < refcoords.size(); i++)
            {
                arc[t][iarc[t]+i] = refcoords[i];
                if (withphysicals)
                    apc[t][iarc[t]+i] = physcoords[i];
            }
            
            iarc[t] += refcoords.size();
        }
        
        if (ln == numleaves-1)
            break;
            
        next();
    }
}

void htracker::getadaptedcoordinates(std::vector<std::vector<double>>& ac, std::vector<std::vector<int>>& leafnums, std::vector<double> noisethreshold)
{
    std::vector<int> ne(8);
    std::vector<std::vector<double>> curvedcoords(8);
    for (int i = 0; i < 8; i++)
    {
        ne[i] = mycurvedelems[i].countedges();
        lagrangeformfunction lff(i,originalcurvatureorder,{});
        curvedcoords[i] = lff.getnodecoordinates();
    }


    // Get the reference ('arc') and physical ('ac') element corner coordinates after all fullsplit adaptation:
    std::vector<std::vector<double>> cornerac;
    std::vector<std::vector<double>> cornerarc;
    atleaves(cornerarc, cornerac, true);
    
    
    // Assign unique edge numbers and deduce edge splits:
    std::vector<int> edgenumbers;
    std:vector<bool> isedgesplit;
    myalgorithm::assignedgenumbers(cornerac, edgenumbers, isedgesplit, noisethreshold);
    

    // Preallocate output containers:
    ac = std::vector<std::vector<double>>(8, std::vector<double>(0));
    transitionsrefcoords = std::vector<std::vector<double>>(8, std::vector<double>(0));
    leavesoftransitions = std::vector<std::vector<int>>(8, std::vector<int>(0));
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
        leavesoftransitions[i] = std::vector<int>(ac[i].size()/ncn[i]/3);
        transitionsrefcoords[i] = std::vector<double>(nn[i]*ac[i].size()/ncn[i]);
    }
    
    resetcursor();
    
    std::vector<double> oc;
    
    int ln = -1; // leaf number
    std::vector<int> iarc(8,0); // indexes in arc
    std::vector<int> firstedge(8,0); // first edge in working element
    for (int i = 0; i < 7; i++)
        firstedge[i+1] = firstedge[i] + ne[i] * cornerarc[i].size()/nn[i]/3;
    std::vector<int> iac(8,0); // index in ac
    while (true)
    {
        if (currentdepth == 0 && indexesinclusters[0] == 0)
            oc = myoriginalelements->getnodecoordinates(parenttypes[0], origindexintype);
    
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
        std::vector<std::vector<int>> splitrefnums = myelems[t].split(splitnum, curedgenums);
      
        // Loop on all transition elements:
        for (int si = 0; si < 8; si++)
        {
            if (splitrefnums[si].size() == 0)
                continue;
        
            std::vector<double> splitrefcoords;
            myelems[t].numstorefcoords(splitrefnums[si], splitrefcoords);
        
            for (int se = 0; se < splitrefcoords.size()/nn[si]/3; se++)
            {
                // Get the ref. coords. of the current transition element:
                std::vector<double> curcoords(3*nn[si]);
                for (int i = 0; i < 3*nn[si]; i++)
                    curcoords[i] = splitrefcoords[se*nn[si]*3+i];
        
                // Bring inside the untransitioned element (if split at all):
                curcoords = myelems[t].calculatecoordinates(curcoords, cornerarc[t], iarc[t], splitnum == 0);
                
                for (int i = 0; i < curcoords.size(); i++)
                    transitionsrefcoords[si][nn[si]*iac[si]/ncn[si]+i] = curcoords[i];
                
                // Make curved:
                if (originalcurvatureorder > 1)
                    curcoords = myelems[si].calculatecoordinates(curvedcoords[si], curcoords);
                    
                // Calculate actual coordinates: 
                curcoords = mycurvedelems[ot].calculatecoordinates(curcoords, oc, 0);

                for (int i = 0; i < curcoords.size(); i++)
                    ac[si][iac[si]+i] = curcoords[i];
                    
                leavesoftransitions[si][iac[si]/ncn[si]/3] = ln;
                
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
            transitionsrefcoords[i].resize(nn[i]*iac[i]/ncn[i]);
            leavesoftransitions[i].resize(iac[i]/ncn[i]/3);
        }
        // ELSE RESIZE TO 0 OR ALREADY EMPTY?
    }
    leafnums = leavesoftransitions;
}

void htracker::tooriginal(std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& oad, std::vector<double>& orc)
{
    std::vector<int> oen;
    getoriginalelementnumber(oen);

    // Preallocate:
    int numrc = 0, no = 0;
    for (int i = 0; i < 8; i++)
    {
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
            int origelem = oen[leavesoftransitions[i][j]];
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
                
            currefs = myelems[i].calculatecoordinates(currefs, transitionsrefcoords[i], 3*nn[i]*j);
    
            int origelem = oen[leavesoftransitions[i][j]];
            for (int k = 0; k < 3*nr; k++)
                orc[oad[origelem]+index[origelem]+k] = currefs[k];
       
            index[origelem] += 3*nr;
        }
    }
}

void htracker::fromoriginal(std::vector<int>& oad, std::vector<double>& orc, std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc)
{
    // Find the transition element in which each ref. coord is:
    int numtran = 0;
    std::vector<int> transindex(8,0);
    for (int i = 0; i < 8; i++)
    {
        transindex[i] = numtran;
        numtran += leavesoftransitions[i].size();
    }    
    // In which transition element is each ref coord?
    std::vector<int> transitionelemnum(orc.size()/3);
    
    // Indexes of the orc's that are in the current tree position:
    std::vector<std::vector<std::vector<int>>> activercs(maxdepth+1, std::vector<std::vector<int>>(10));
    
    resetcursor(true);
    
    int origelem = -1;
    int ln = -1; // leaf number
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ot = parenttypes[0];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        if (ns == 0)
        {
            origelem++;
            int numrefsinorig = (oad[origelem+1]-oad[origelem])/3;
            activercs[0] = {std::vector<int>(numrefsinorig)};
            for (int i = 0; i < numrefsinorig; i++)
                activercs[0][0][i] = oad[origelem]/3+i;
        }
    
        if (isatleaf())
        {
            ln++;
         
        }
        else
        {
            if (ns > 0)
            {
                std::vector<int> par = activercs[ns-1][indexesinclusters[currentdepth-1]];
                int numinparent = par.size();
                
                std::vector<double> parcoords(3*numinparent);
                for (int i = 0; i < numinparent; i++)
                {
                    parcoords[3*i+0] = orc[3*par[i]+0];
                    parcoords[3*i+1] = orc[3*par[i]+1];
                    parcoords[3*i+2] = orc[3*par[i]+2];
                }
            
                std::vector<double> refcoords = getreferencecoordinates();
                
                // Redirect the ref. coord. to the current element if inside it:
                std::vector<bool> isinside;
                myelems[t].isinsideelement(parcoords, refcoords, isinside, 1e-10);
                
                myalgorithm::splitvector(par, isinside, activercs[ns-1][indexesinclusters[currentdepth-1]], activercs[ns][ic]);
            }
        }
        
        if (ln == numleaves-1)
            break;
            
        next();
    }
}

void htracker::tostorage(void)
{
    transitionsrefcoords = {};
}

void htracker::print(void)
{
    std::cout << "#leaves is " << numleaves << " | max #splits is " << maxdepth << std::endl;

    for (int i = 0; i < splitdata.size(); i++)
        std::cout << splitdata[i];
        
    std::cout << " (" << splitdata.size() << " bits)" << std::endl;
}

int htracker::countbits(void)
{
    return splitdata.size();
}

