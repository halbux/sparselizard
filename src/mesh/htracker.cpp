#include "htracker.h"
#include "myalgorithm.h"
#include "lagrangeformfunction.h"


htracker::htracker(std::shared_ptr<rawmesh> origmesh, int curvatureorder, std::vector<int> numelemspertype)
{
    myoriginalmesh = origmesh;

    originalcurvatureorder = curvatureorder;
    originalcount = numelemspertype;
    
    if (curvatureorder == -1)
        originalcurvatureorder = origmesh->getelements()->getcurvatureorder();
    if (numelemspertype.size() == 0)
    {
        int dim = origmesh->getelements()->getdimension();
        originalcount = std::vector<int>(8,0);
        for (int i = 0; i < 8; i++)
        {
            element myelem(i);
            if (myelem.getelementdimension() == dim)
                originalcount[i] = origmesh->getelements()->count(i);
        }
    }
        
    myelems = std::vector<element>(8);
    mycurvedelems = std::vector<element>(8);
    nn = std::vector<int>(8);
    ncn = std::vector<int>(8);
    straightrefcoords = std::vector<std::vector<double>>(8);
    curvedrefcoords = std::vector<std::vector<double>>(8);
    for (int i = 0; i < 8; i++)
    {
        numleaves += originalcount[i];
        lagrangeformfunction lff(i,1,{});
        straightrefcoords[i] = lff.getnodecoordinates();
        // LFF NOT YET DEFINED FOR CURVED PYRAMIDS
        if (i < 7)
        {
            lagrangeformfunction lffc(i,originalcurvatureorder,{});
            curvedrefcoords[i] = lffc.getnodecoordinates();
        }
        myelems[i] = element(i);
        mycurvedelems[i] = element(i,originalcurvatureorder);
        nn[i] = myelems[i].countnodes();
        ncn[i] = mycurvedelems[i].countcurvednodes();
    }
    
    // Initialize the tree:
    splitdata = std::vector<bool>(numleaves, false);
    
    // Initialize the remaining containers:
    transitionsrefcoords = std::vector<std::vector<double>>(8, std::vector<double>(0));
    leavesoftransitions = std::vector<std::vector<int>>(8, std::vector<int>(0));
    originalsoftransitions = std::vector<std::vector<int>>(8, std::vector<int>(0));
    touser = std::vector<std::vector<int>>(8, std::vector<int>(0));
    toht = std::vector<std::vector<int>>(8, std::vector<int>(0));

    int ind = 0;
    for (int i = 0; i < 8; i++)
    {
        int num = originalcount[i];
        transitionsrefcoords[i] = myalgorithm::duplicate(straightrefcoords[i], num);
        // Leaves equal originals here:
        leavesoftransitions[i] = myalgorithm::getequallyspaced(ind, 1, num);
        originalsoftransitions[i] = std::vector<int>(2*num);
        for (int j = 0; j < num; j++)
        {
            originalsoftransitions[i][2*j+0] = i;
            originalsoftransitions[i][2*j+1] = j;
        }
        touser[i] = myalgorithm::getequallyspaced(0, 1, num);
        toht[i] = myalgorithm::getequallyspaced(0, 1, num);
        
        ind += num;
    }
}

std::shared_ptr<rawmesh> htracker::getoriginalmesh(void)
{
    if (myoriginalmesh.expired())
    {
        std::cout << "Error in 'htracker' object: the original mesh is needed but it was destroyed" << std::endl;
        abort();
    }
    return myoriginalmesh.lock();
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
    if (isrefcalc)
        parentrefcoords[0] = {straightrefcoords[parenttypes[0]]};
}

int htracker::next(void)
{
    elements* myelements = getoriginalmesh()->getelements();

    int throughedgenum = -1;

    // If not split we are at a leaf:
    if (splitdata[cursorposition] == false)
    {
        cursorposition++;
        
        if (currentdepth > 0)
        {            
            // Move back up in the tree as much as needed:
            indexesinclusters[currentdepth]++;
            while (currentdepth > 0 && indexesinclusters[currentdepth] == numsubelems[parenttypes[currentdepth-1]])
            {
                currentdepth--;
                indexesinclusters[currentdepth]++;
            }
            // For pyramids the subelement type might change from 4 to 7:
            if (currentdepth > 0 && parenttypes[currentdepth-1] == 7 && indexesinclusters[currentdepth] == 4)
                parenttypes[currentdepth] = 7;
        }
        if (currentdepth == 0) // not else here because currentdepth might change above
        {
            indexesinclusters[0] = 0;
            
            origindexintype++;
            // Skip empty types:
            while (origindexintype == originalcount[parenttypes[0]])
            {
                parenttypes[0]++;
                origindexintype = 0;
            }
            if (isrefcalc)
                parentrefcoords[0] = {straightrefcoords[parenttypes[0]]};
        }
    }
    else
    {
        cursorposition++;
        
        // Get the through-edge number (if any):
        int t = parenttypes[currentdepth];
        int ic = indexesinclusters[currentdepth];
        if (t == 4)
        {
            throughedgenum = 0;
            if (splitdata[cursorposition])
                throughedgenum += 2;
            if (splitdata[cursorposition+1])
                throughedgenum += 1;
                
            // If it is 3 we have to define it:
            if (isrefcalc && throughedgenum == 3)
            {
                // Get the physical node coordinates of the original element:
                std::vector<double> origelemcoords = myelements->getnodecoordinates(parenttypes[0], origindexintype);
                // Get the curved reference coordinates in the original element:
                std::vector<double> refsinorig;
                if (currentdepth > 0)
                    refsinorig = myelems[t].calculatecoordinates(curvedrefcoords[t], parentrefcoords[currentdepth][ic], 0, originalcurvatureorder == 1);
                // Calculate the physical coordinates of the current element:
                std::vector<double> coordsinorigelem = mycurvedelems[parenttypes[0]].calculatecoordinates(refsinorig, origelemcoords, 0, currentdepth == 0);
                // Calculate the best through-edge number:
                throughedgenum = mycurvedelems[4].choosethroughedge(coordsinorigelem);
                
                // Write it to the tree (status is 11 now):
                if (throughedgenum%2 == 0)
                    splitdata[cursorposition+1] = false;
                if (((throughedgenum-throughedgenum%2)/2)%2 == 0)
                    splitdata[cursorposition] = false;
            }
            cursorposition += 2;
        }
        
        // Optionally compute the subelements reference coordinates in the original element:
        if (isrefcalc)
        {
            std::vector<std::vector<double>> cornerrefcoords;
            myelems[t].fullsplit(cornerrefcoords, throughedgenum);
            
            int ind = 0;
            for (int i = 0; i < 8; i++)
            {
                for (int j = 0; j < cornerrefcoords[i].size()/nn[i]/3; j++)
                {
                    std::vector<double> cc(3*nn[i]);
                    for (int k = 0; k < cc.size(); k++)
                        cc[k] = cornerrefcoords[i][j*3*nn[i]+k];
                    
                    parentrefcoords[currentdepth+1][ind] = myelems[t].calculatecoordinates(cc, parentrefcoords[currentdepth][ic]);

                    ind++;
                }
            }
        }
    
        // Split pyramids give 4 tets and 6 pyramids:
        if (parenttypes[currentdepth] != 7)
            parenttypes[currentdepth+1] = parenttypes[currentdepth];
        else
            parenttypes[currentdepth+1] = 4;

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

std::vector<double> htracker::getreferencecoordinates(void)
{
    return parentrefcoords[currentdepth][indexesinclusters[currentdepth]];
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

void htracker::getoriginalelement(std::vector<int>& oet, std::vector<int>& oei)
{
    oet = std::vector<int>(numleaves);
    oei = std::vector<int>(numleaves);
    
    resetcursor();
    
    int ln = -1;
    while (true)
    { 
        if (isatleaf())
        {
            ln++;
            oet[ln] = parenttypes[0];
            oei[ln] = origindexintype;
        }
        
        if (ln == numleaves-1)
            break;
        
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

int htracker::counttransitions(int elementtypenumber)
{
    return leavesoftransitions[elementtypenumber].size();
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

    while (true)
    {
        int t = parenttypes[currentdepth];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
        
        // Reset cluster info:
        if (ic == 0)
            ngr[ns] = 0;

        if (not(isatleaf()))
        {
            newsplitdata[ni] = true;
            ni++;
            
            // Write the through-edge number (if any):
            if (t == 4)
            {
                newsplitdata[ni] = splitdata[cursorposition+1];
                newsplitdata[ni+1] = splitdata[cursorposition+2];
                ni += 2;
            }
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
                // Initialize through-edge number at 3 (undefined):
                if (t == 4)
                {
                    newsplitdata[ni] = true;
                    newsplitdata[ni+1] = true;
                    ni += 2;
                }
                for (int i = 0; i < numsubelems[t]; i++)
                    newsplitdata[ni+i] = false;
                ni += numsubelems[t];
                
                if (newmaxdepth < ns+1)
                    newmaxdepth = ns+1;
            }
            else
            {
                newsplitdata[ni] = false;
                ni++;
        
                if (operations[ln] == -1)
                    ngr[ns]++;
                
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
        }
        
        if (ln == numleaves-1)
            break;
            
        next();
    }

    newsplitdata.resize(ni);
    splitdata = newsplitdata;
    numleaves = newnumleaves;
    maxdepth = newmaxdepth;
}

void htracker::atleaves(std::vector<std::vector<double>>& arc, std::vector<std::vector<double>>& apc, bool withphysicals)
{
    elements* myelements = getoriginalmesh()->getelements();

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
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        if (withphysicals && ns == 0)
            oc = myelements->getnodecoordinates(t, origindexintype);
    
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

void htracker::getadaptedcoordinates(std::vector<std::vector<double>>& ac)
{
    elements* myelements = getoriginalmesh()->getelements();
    
    std::vector<int> ne(8);
    for (int i = 0; i < 8; i++)
        ne[i] = myelems[i].countedges();

    // Get the reference ('arc') and physical ('apc') element corner coordinates after all fullsplit adaptations:
    std::vector<std::vector<double>> cornerarc, cornerapc;
    atleaves(cornerarc, cornerapc, true);    
    
    // Assign unique edge numbers and deduce edge splits:
    std::vector<int> edgenumbers;
    std:vector<bool> isedgesplit;
    myalgorithm::assignedgenumbers(cornerapc, edgenumbers, isedgesplit);
    

    // Preallocate output containers to upper bound size:
    ac = std::vector<std::vector<double>>(8, std::vector<double>(0));
    transitionsrefcoords = std::vector<std::vector<double>>(8, std::vector<double>(0));
    leavesoftransitions = std::vector<std::vector<int>>(8, std::vector<int>(0));
    originalsoftransitions = std::vector<std::vector<int>>(8, std::vector<int>(0));

    std::vector<int> ub = countupperbound();
    for (int i = 0; i < 8; i++)
    {
        ac[i] = std::vector<double>(3*ncn[i]*ub[i]);
        leavesoftransitions[i] = std::vector<int>(ub[i]);
        originalsoftransitions[i] = std::vector<int>(2*ub[i]);
        transitionsrefcoords[i] = std::vector<double>(3*nn[i]*ub[i]);
    }
    
    resetcursor();
    
    std::vector<double> oc;
    
    int ln = -1; // leaf number
    std::vector<int> iarc(8,0); // indexes in arc
    std::vector<int> firstedge(8,0); // first edge in working element
    for (int i = 0; i < 7; i++)
        firstedge[i+1] = firstedge[i] + ne[i] * cornerarc[i].size()/nn[i]/3;
    std::vector<int> ite(8,0); // index of trans. elem.
    while (true)
    {
        int t = parenttypes[currentdepth];

        if (currentdepth == 0)
            oc = myelements->getnodecoordinates(t, origindexintype);
    
        while (not(isatleaf()))
            next();
    
        ln++;
        
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
                    transitionsrefcoords[si][3*nn[si]*ite[si]+i] = curcoords[i];
                
                // Make curved:
                if (originalcurvatureorder > 1)
                    curcoords = myelems[si].calculatecoordinates(curvedrefcoords[si], curcoords);
                    
                // Calculate actual coordinates: 
                curcoords = mycurvedelems[parenttypes[0]].calculatecoordinates(curcoords, oc, 0);

                for (int i = 0; i < curcoords.size(); i++)
                    ac[si][3*ncn[si]*ite[si]+i] = curcoords[i];
                    
                leavesoftransitions[si][ite[si]] = ln;
                originalsoftransitions[si][2*ite[si]+0] = parenttypes[0];
                originalsoftransitions[si][2*ite[si]+1] = origindexintype;
                
                ite[si]++;
                
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
        ac[i].resize(3*ncn[i]*ite[i]);
        transitionsrefcoords[i].resize(3*nn[i]*ite[i]);
        leavesoftransitions[i].resize(ite[i]);
        originalsoftransitions[i].resize(2*ite[i]);
        
        touser[i] = myalgorithm::getequallyspaced(0, 1, ite[i]);
        toht[i] = myalgorithm::getequallyspaced(0, 1, ite[i]);
    }
}

std::vector<int> htracker::countupperbound(void)
{
    std::vector<int> nit = countintypes();
    
    std::vector<int> output(8,0);
    
    output[0] = nit[0];
    output[1] = 2*nit[1];
    output[2] = 4*nit[2] + 4*nit[3];
    output[3] = 4*nit[3];
    output[4] = 8*nit[4];
    // DEFINE ALSO FOR HEXAHEDRA, PRISMS AND PYRAMIDS (+ adapt tet nums)
    
    return output;
}

void htracker::renumbertransitions(std::vector<std::vector<int>>& renumbering)
{
    for (int i = 0; i < 8; i++)
    {
        touser[i] = myalgorithm::chainrenumbering(touser[i], renumbering[i]);
        toht[i] = myalgorithm::invertrenumbering(touser[i]);
    }
}

void htracker::tooriginal(std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& oad, std::vector<double>& orc, std::vector<std::vector<int>>& maprctoorc)
{
    std::vector<int> oen;
    getoriginalelementnumber(oen);

    // Preallocate:
    maprctoorc = std::vector<std::vector<int>>(8, std::vector<int>(0));
    for (int i = 0; i < 8; i++)
        maprctoorc[i] = std::vector<int>(rc[i].size()/3);
        
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
            for (int k = 0; k < nr; k++)
                maprctoorc[i][ad[i][j]/3+k] = (oad[origelem]+index[origelem])/3+k;
       
            index[origelem] += 3*nr;
        }
    }
}

void htracker::fromoriginal(std::vector<int>& oad, std::vector<double>& orc, std::vector<std::vector<int>>& ad, std::vector<std::vector<double>>& rc, std::vector<int>& maporctorc)
{
    int elemdim = getoriginalmesh()->getelements()->getdimension();

    ad = std::vector<std::vector<int>>(8, std::vector<int>(0));
    rc = std::vector<std::vector<double>>(8, std::vector<double>(0));
    for (int i = 0; i < 8; i++)
    {
        int curnumtran = leavesoftransitions[i].size();
        if (curnumtran > 0)
        {
            ad[i] = std::vector<int>(curnumtran+1,0);
            // Max bound (size is unknown for now):
            rc[i] = std::vector<double>(orc.size());
        }
    }    
    maporctorc = std::vector<int>(2*orc.size()/3, -1);
    
    // Indexes of the 'orc' points that are in the current tree position:
    std::vector<std::vector<int>> actives(maxdepth+1, std::vector<int>(0));
    
    // Needed for the root finding:
    std::vector<polynomials> polys(8);
    for (int i = 0; i < 8; i++)
        polys[i] = polynomials(lagrangeformfunction(i,1,{}).getformfunctionpolynomials());
    
    resetcursor(true);
    
    int origelem = -1;
    int ln = -1; // leaf number
    std::vector<int> ti(8,0); // transition element index
    
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        // Update 'actives':
        if (ns == 0)
        {
            origelem++;
            // Set all to active:
            int numrefsinorig = (oad[origelem+1]-oad[origelem])/3;
            actives[0] = myalgorithm::getequallyspaced(oad[origelem]/3, 1, numrefsinorig);
        }
        else
        {
            std::vector<double> refcoords = getreferencecoordinates();
        
            // Actives in parent:
            std::vector<int> par = actives[ns-1];
            int numactivesinparent = par.size();
            
            std::vector<double> parcoords(3*numactivesinparent);
            for (int i = 0; i < numactivesinparent; i++)
            {
                parcoords[3*i+0] = orc[3*par[i]+0];
                parcoords[3*i+1] = orc[3*par[i]+1];
                parcoords[3*i+2] = orc[3*par[i]+2];
            }
            
            // Redirect the reference coordinates to the current element if inside it:
            std::vector<bool> isinside;
            myelems[t].isinsideelement(parcoords, refcoords, isinside, 1e-10);
            
            myalgorithm::splitvector(par, isinside, actives[ns-1], actives[ns]);
        }
    
        if (isatleaf())
        {
            ln++;
            
            // Loop on all transition elements of this leaf:
            for (int i = 0; i < 8; i++)
            {
                while (ti[i] < leavesoftransitions[i].size() && leavesoftransitions[i][ti[i]] == ln)
                {
                    // Get the reference coordinates of the current transition element:
                    std::vector<double> currefcoords(3*nn[i]);
                    for (int j = 0; j < 3*nn[i]; j++)
                        currefcoords[j] = transitionsrefcoords[i][3*nn[i]*ti[i]+j];
                        
                    // Actives in the current leaf:
                    std::vector<int> activesinleaf = actives[ns];
                    int numactivesinleaf = activesinleaf.size();
                    
                    std::vector<double> activecoords(3*numactivesinleaf);
                    for (int j = 0; j < numactivesinleaf; j++)
                    {
                        activecoords[3*j+0] = orc[3*activesinleaf[j]+0];
                        activecoords[3*j+1] = orc[3*activesinleaf[j]+1];
                        activecoords[3*j+2] = orc[3*activesinleaf[j]+2];
                    }
                    
                    std::vector<bool> isintrans;
                    myelems[i].isinsideelement(activecoords, currefcoords, isintrans, 1e-10);
                    // Actives in the current transition element:
                    std::vector<int> activesintrans;
                    myalgorithm::splitvector(activesinleaf, isintrans, actives[ns], activesintrans);
            
                    // Populate 'ad':
                    ad[i][ti[i]+1] = ad[i][ti[i]]+3*activesintrans.size();
            
                    if (activesintrans.size() > 0)
                    {
                        // Find the corresponding reference coordinate in the transition element's own reference.
                        // First create the polynomials for the system to solve.
                        std::vector<double> xyz = myalgorithm::separate(currefcoords, 3, myalgorithm::getequallyspaced(0,1,elemdim));
                        polynomials syspolys = polys[i].sum(xyz);
                
                        // Loop on all actives:
                        for (int j = 0; j < activesintrans.size(); j++)
                        {                        
                            std::vector<double> kietaphi = {0.0,0.0,0.0};
                            std::vector<double> rhs = {orc[3*activesintrans[j]+0], orc[3*activesintrans[j]+1], orc[3*activesintrans[j]+2]};
                            
                            if (myalgorithm::getroot(syspolys, rhs, kietaphi) == 1 && myelems[i].isinsideelement(kietaphi[0], kietaphi[1], kietaphi[2]))
                            {
                                rc[i][ad[i][ti[i]]+3*j+0] = kietaphi[0];
                                rc[i][ad[i][ti[i]]+3*j+1] = kietaphi[1];
                                rc[i][ad[i][ti[i]]+3*j+2] = kietaphi[2];
                                
                                maporctorc[2*activesintrans[j]+0] = i;
                                maporctorc[2*activesintrans[j]+1] = ad[i][ti[i]]/3+j;
                            }
                            else
                            {
                                std::cout << "Error in 'htracker' object: root finding algorithm for mesh adaptivity failed to converge for a " << myelems[i].gettypename() << " element" << std::endl;
                                abort();
                            }
                        }
                    }
                    
                    ti[i]++;
                }
            }
        }
        
        if (ln == numleaves-1)
            break;
            
        next();
    }
    
    // Fit to actual size:
    for (int i = 0; i < 8; i++)
    {
        if (ad[i].size() > 0)
            rc[i].resize(ad[i][ad[i].size()-1]);
    }
    
    // Sanity check:
    for (int i = 0; i < maporctorc.size(); i++)
    {
        if (maporctorc[i] < 0)
        {
            std::cout << "Error in 'htracker' object: mapping failed for at least one point" << std::endl;
            abort();
        }
    }
}

void htracker::getattarget(std::vector<std::vector<int>>& userad, std::vector<std::vector<double>>& userrc, htracker* target, std::vector<std::vector<int>>& targettranselems, std::vector<std::vector<double>>& targetrefcoords)
{
    // Take the transition element renumbering into account:
    std::vector<std::vector<int>> ad(8, std::vector<int>(0));
    std::vector<std::vector<double>> rc(8, std::vector<double>(0));
    for (int i = 0; i < 8; i++)
    {
        if (userrc[i].size() > 0)
            myalgorithm::reorder(userad[i], userrc[i], toht[i], ad[i], rc[i]);
    }


    // Get from this htracker's transition elements to the original elements:
    std::vector<int> oad;
    std::vector<double> orc;
    std::vector<std::vector<int>> maprctoorc;
    
    tooriginal(ad, rc, oad, orc, maprctoorc);
    
    // Get from the original elements to the target htracker's transition elements:
    std::vector<std::vector<int>> tad;
    std::vector<std::vector<double>> trc;
    std::vector<int> maporctorc;
    
    target->fromoriginal(oad, orc, tad, trc, maporctorc);
    
    
    // Create output containers:
    targettranselems = std::vector<std::vector<int>>(8, std::vector<int>(0));
    targetrefcoords = std::vector<std::vector<double>>(8, std::vector<double>(0));
    
    // Create a vector to map the position in the target 'rc' to the target transition element number:
    std::vector<std::vector<int>> map(8, std::vector<int>(0));
    for (int i = 0; i < 8; i++)
    {
        if (trc[i].size() == 0)
            continue;
    
        map[i] = std::vector<int>(trc[i].size()/3);
    
        int index = 0;
        for (int j = 0; j < tad[i].size()-1; j++)
        {
            for (int k = 0; k < (tad[i][j+1]-tad[i][j])/3; k++)
            {
                map[i][index] = j;
                index++;
            }
        }
    }
    
    // Loop on all 'rc' entries:
    for (int i = 0; i < 8; i++)
    {
        int nr = rc[i].size()/3;
        
        if (nr == 0)
            continue;
        
        targettranselems[i] = std::vector<int>(2*nr);
        targetrefcoords[i] = std::vector<double>(3*nr);
        
        for (int j = 0; j < ad[i].size()-1; j++)
        {
            int pos = ad[i][j]/3;
            int userpos = userad[i][touser[i][j]]/3;
            for (int k = 0; k < (ad[i][j+1]-ad[i][j])/3; k++)
            {
                // To original element:
                int orc = maprctoorc[i][pos+k];
                
                // To target:
                int ttype = maporctorc[2*orc+0];
                int tind = maporctorc[2*orc+1];
            
                targettranselems[i][2*(userpos+k)+0] = ttype;
                targettranselems[i][2*(userpos+k)+1] = target->touser[ttype][map[ttype][tind]];
                
                targetrefcoords[i][3*(userpos+k)+0] = trc[ttype][3*tind+0];
                targetrefcoords[i][3*(userpos+k)+1] = trc[ttype][3*tind+1];
                targetrefcoords[i][3*(userpos+k)+2] = trc[ttype][3*tind+2];
            }
        }
    }
}

void htracker::atoriginal(int tt, int tn, int& ort, int& orn, std::vector<int>& on, std::vector<int>& oe, std::vector<int>& of)
{
    tn = toht[tt][tn];

    ort = originalsoftransitions[tt][2*tn+0];
    orn = originalsoftransitions[tt][2*tn+1];

    std::vector<double> nc(3*nn[tt]);
    for (int i = 0; i < 3*nn[tt]; i++)
        nc[i] = transitionsrefcoords[tt][3*nn[tt]*tn+i];
    myelems[ort].atnode(nc, on);
    
    if (tt > 1)
    {
        std::vector<double> eb = myelems[tt].getedgebarycenter(nc);
        myelems[ort].atedge(eb, oe);
    }
    if (tt > 3)
    {
        std::vector<double> fb = myelems[tt].getfacebarycenter(nc);
        myelems[ort].atface(fb, of);
    }
}

void htracker::getattarget(std::vector<int>& olv, htracker* target, std::vector<int>& tlv)
{
    tlv = std::vector<int>(target->numleaves, -1);

    resetcursor();
    target->resetcursor();
    
    int oln = 0;
    int tln = 0;
    
    while (true)
    {
        // Move to next leaf:
        while (not(isatleaf()))
            next();
        while (not(target->isatleaf()))
            target->next();
            
        int ocd = currentdepth;
        int tcd = target->currentdepth;
        
        int oic = indexesinclusters[ocd];
        int tic = target->indexesinclusters[tcd];
        
        int ons = -1;
        int tns = -1;
        if (ocd > 0)
            ons = numsubelems[parenttypes[ocd-1]];
        if (tcd > 0)
            tns = numsubelems[target->parenttypes[tcd-1]];
        
        // If both leaves match:
        if (ocd == tcd)
            tlv[tln] = olv[oln];
        // If original tree is one deeper here:
        if (ocd > tcd)
            tlv[tln] = std::max(olv[oln], tlv[tln]);
        // If target tree is one deeper here:
        if (ocd < tcd)
            tlv[tln] = olv[oln];
        
        if (ocd >= tcd || tic == tns-1)
        {
            // Will be reached for original and target at the same time:
            if (oln == numleaves-1)
                break;
        
            next();
            oln++;
        }
        if (ocd <= tcd || oic == ons-1)
        {
            target->next();
            tln++;
        }
    }
}

int htracker::getleafnumber(int transitiontype, int transitionnumber)
{
    int tn = toht[transitiontype][transitionnumber];
    return leavesoftransitions[transitiontype][tn];
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

