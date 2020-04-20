#include "htracker.h"
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
        else
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

    // Is any element split/tagged for splitting in the cluster?
    std::vector<bool> isanysplit(maxdepth+1);
    // Is any element tagged for grouping in the cluster?
    std::vector<bool> isanygroup(maxdepth+1);
    
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
        {
            isanysplit[ns] = false;
            isanygroup[ns] = false;
        }

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
            isanysplit[ns] = true;
            newsplitdata[ni] = true;
            ni++;
        }
        else
        {
            ln++;
            
            // Split:
            if (operations[ln] == 1)
            {
                isanysplit[ns] = true;
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
                
                if (newmaxdepth <= ns)
                    newmaxdepth = ns+1;
            }
            
            if (operations[ln] == -1)
                isanygroup[ns] = true;
                
            // Group/unchanged:
            if (operations[ln] < 1)
            {
                newsplitdata[ni] = false;
                ni++;
            }
            
            // A cluster to group must always end with a leaf:
            if (ns > 0 && ic == numsubelems[parenttypes[ns-1]]-1)
            {
                // Group the cluster (if allowed):
                int md = ns;
                if (not(isanysplit[ns]) && isanygroup[ns])
                {
                    int pt = parenttypes[ns-1];
                    newnumleaves -= numsubelems[pt]-1;
                    ni -= numsubelems[pt];
                    if (pt == 4)
                        ni -= 2;
                    newsplitdata[ni-1] = false;
                    md--;
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

void htracker::getadaptedrefcoords(std::vector<std::vector<double>>& arc)
{
    std::vector<int> nit = countintypes();
    std::vector<int> nn = {1,2,3,4,4,8,6,5}; // number of corner nodes

    // Preallocate:
    std::vector<element> els(8);
    arc = std::vector<std::vector<double>>(8, std::vector<double>(0));
    std::vector<std::vector<double>> crc(8); // corner ref coords
    for (int i = 0; i < 8; i++)
    {
        arc[i] = std::vector<double>(3*nn[i]*nit[i]);
        els[i] = element(i);
        lagrangeformfunction lff(i,1,{});
        crc[i] = lff.getnodecoordinates();
    }

    // parentrc[depth][indexincluster]:
    std::vector<std::vector<std::vector<double>>> parentrc(maxdepth+1, std::vector<std::vector<double>>(10));
    
    resetcursor();
    
    int ln = -1;
    std::vector<int> iarc(8,0); // indexes in arc
    while (true)
    {
        int t = parenttypes[currentdepth];
        int ns = currentdepth;
        int ic = indexesinclusters[currentdepth];
    
        if (ns == 0)
            parentrc[0] = {crc[t]};
    
        if (isatleaf())
        {
            ln++;
            
            for (int i = 0; i < parentrc[ns][ic].size(); i++)
                arc[t][iarc[t]+i] = parentrc[ns][ic][i];
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

