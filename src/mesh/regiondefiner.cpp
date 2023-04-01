#include "regiondefiner.h"


void regiondefiner::defineskinregion(int regnum)
{
    bool isnotall = (toskin[regnum] != -1);

    if (isnotall)
        myphysicalregions->errorundefined({toskin[regnum]});

    physicalregion* newphysreg = myphysicalregions->get(skins[regnum]);
    physicalregion* curphysreg;
    if (isnotall)
        curphysreg = myphysicalregions->get(toskin[regnum]);

    int physregdim;
    if (isnotall)
        physregdim = curphysreg->getelementdimension();
    else
        physregdim = myelements->getdimension();
        
    if (physregdim == 0)
    {
        logs log;
        log.msg() << "Error in 'regiondefiner' object: cannot get the skin of point elements" << std::endl;
        log.error();
    }

    // The skin has a one lower dimension:
    int skindim = physregdim-1;

    std::vector<std::vector<int>>* curelems;
    if (isnotall)
        curelems = curphysreg->getelementlist();

    // Loop on all skin element types:
    for (int skinelemtype = 0; skinelemtype <= 7; skinelemtype++)
    {
        element myskinelement(skinelemtype);
        if (myskinelement.getelementdimension() != skindim)
            continue;

        // Vector to count the number of times a skin element appears in an element:
        std::vector<int> numoccurences(myelements->count(skinelemtype),0);

        // Loop on all element types in the physical region:
        for (int elemtype = 0; elemtype <= 7; elemtype++)
        {
            int numelems;
            if (isnotall)
                numelems = curelems->at(elemtype).size();
            else
                numelems = myelements->countcells(elemtype);
            
            if (numelems == 0)
                continue;

            element myelement(elemtype);

            int numsubtypeelems = myelement.counttype(skinelemtype);
            if (numsubtypeelems == 0)
                continue;

            // Loop on all elements:
            for (int elem = 0; elem < numelems; elem++)
            {
                int curelem;
                if (isnotall)
                    curelem = curelems->at(elemtype)[elem];
                else
                    curelem = elem;

                for (int subelem = 0; subelem < numsubtypeelems; subelem++)
                    numoccurences[myelements->getsubelement(skinelemtype, elemtype, curelem, subelem)]++;
            }
        }

        // All single occurences are skin elements:
        for (int j = 0; j < numoccurences.size(); j++)
        {
            if (numoccurences[j] == 1)
                newphysreg->addelement(skinelemtype, j);
        }
    }
    newphysreg->removeduplicatedelements();
}

void regiondefiner::defineboxregion(int regnum)
{
    bool isnotall = (tobox[regnum] != -1);
    
    if (isnotall)
        myphysicalregions->errorundefined({tobox[regnum]});

    std::vector<double>* nodecoords = mynodes->getcoordinates();

    std::vector<double> boxlimit = boxlimits[regnum];
    
    // Make the box limit slightly larger to remove the roundoff noise issues:
    boxlimit[0] -= noisethreshold; boxlimit[2] -= noisethreshold; boxlimit[4] -= noisethreshold;
    boxlimit[1] += noisethreshold; boxlimit[3] += noisethreshold; boxlimit[5] += noisethreshold;

    physicalregion* newphysreg = myphysicalregions->get(boxed[regnum]);
    physicalregion* curphysreg;
    if (isnotall)
        curphysreg = myphysicalregions->get(tobox[regnum]);

    std::vector<std::vector<int>>* curelems;
    if (isnotall)
        curelems = curphysreg->getelementlist();

    // Loop on all element types:
    for (int elemtype = 0; elemtype <= 7; elemtype++)
    {
        int numelems;
        if (isnotall)
            numelems = curelems->at(elemtype).size();
        else
            numelems = myelements->countcells(elemtype);
        
        if (numelems == 0)
            continue;

        element myelement(elemtype);

        // Loop on all subelement types:
        for (int subelemtype = 0; subelemtype <= 7; subelemtype++)
        {
            // Skip all subelements that are not of the requested dimension:
            element mysubelement(subelemtype);
            int numnodes = mysubelement.countcurvednodes();
            int numsubtypeelems = myelement.counttype(subelemtype);

            if (mysubelement.getelementdimension() != boxelemdims[regnum] || numsubtypeelems == 0)
                continue;

            // Loop on all elements:
            for (int elem = 0; elem < numelems; elem++)
            {
                int curelem;
                if (isnotall)
                    curelem = curelems->at(elemtype)[elem];
                else
                    curelem = elem;

                for (int subelem = 0; subelem < numsubtypeelems; subelem++)
                {
                    int cursubelem = myelements->getsubelement(subelemtype, elemtype, curelem, subelem);

                    // Check if the coordinates of all nodes in the subelement are within the box limits:
                    bool isinlimits = true;
                    for (int node = 0; node < numnodes; node++)
                    {
                        int curnode = myelements->getsubelement(0, subelemtype, cursubelem, node);

                        double curnodex = nodecoords->at(3*curnode+0);
                        double curnodey = nodecoords->at(3*curnode+1);
                        double curnodez = nodecoords->at(3*curnode+2);

                        if (curnodex < boxlimit[0] || curnodex > boxlimit[1] || curnodey < boxlimit[2] || curnodey > boxlimit[3] || curnodez < boxlimit[4] || curnodez > boxlimit[5])
                        {
                            isinlimits = false;
                            break;
                        }
                    }
                    if (isinlimits)
                        newphysreg->addelement(subelemtype, cursubelem);
                }
            }
        }
    }
    newphysreg->removeduplicatedelements();
}

void regiondefiner::definesphereregion(int regnum)
{
    bool isnotall = (tosphere[regnum] != -1);
    
    if (isnotall)
        myphysicalregions->errorundefined({tosphere[regnum]});

    std::vector<double>* nodecoords = mynodes->getcoordinates();

    std::vector<double> spherecenter = spherecenters[regnum];
    double sphereradius = sphereradii[regnum];
    
    // Make the sphere radius slightly larger to remove the roundoff noise issues:
    sphereradius += noisethreshold;

    physicalregion* newphysreg = myphysicalregions->get(sphered[regnum]);
    physicalregion* curphysreg;
    if (isnotall)
        curphysreg = myphysicalregions->get(tosphere[regnum]);

    std::vector<std::vector<int>>* curelems;
    if (isnotall)
        curelems = curphysreg->getelementlist();

    // Loop on all element types:
    for (int elemtype = 0; elemtype <= 7; elemtype++)
    {
        int numelems;
        if (isnotall)
            numelems = curelems->at(elemtype).size();
        else
            numelems = myelements->countcells(elemtype);
        
        if (numelems == 0)
            continue;

        element myelement(elemtype);

        // Loop on all subelement types:
        for (int subelemtype = 0; subelemtype <= 7; subelemtype++)
        {
            // Skip all subelements that are not of the requested dimension:
            element mysubelement(subelemtype);
            int numnodes = mysubelement.countcurvednodes();
            int numsubtypeelems = myelement.counttype(subelemtype);

            if (mysubelement.getelementdimension() != sphereelemdims[regnum] || numsubtypeelems == 0)
                continue;

            // Loop on all elements:
            for (int elem = 0; elem < numelems; elem++)
            {
                int curelem;
                if (isnotall)
                    curelem = curelems->at(elemtype)[elem];
                else
                    curelem = elem;

                for (int subelem = 0; subelem < numsubtypeelems; subelem++)
                {
                    int cursubelem = myelements->getsubelement(subelemtype, elemtype, curelem, subelem);

                    // Check if the coordinates of all nodes in the subelement are within the sphere limits:
                    bool isinlimits = true;
                    for (int node = 0; node < numnodes; node++)
                    {
                        int curnode = myelements->getsubelement(0, subelemtype, cursubelem, node);

                        double curnodex = nodecoords->at(3*curnode+0);
                        double curnodey = nodecoords->at(3*curnode+1);
                        double curnodez = nodecoords->at(3*curnode+2);

                        if (std::sqrt(std::pow(spherecenter[0]-curnodex,2) + std::pow(spherecenter[1]-curnodey,2) + std::pow(spherecenter[2]-curnodez,2)) > sphereradius)
                        {
                            isinlimits = false;
                            break;
                        }
                    }
                    if (isinlimits)
                        newphysreg->addelement(subelemtype, cursubelem);
                }
            }
        }
    }
    newphysreg->removeduplicatedelements();
}

void regiondefiner::defineexclusionregion(int regnum)
{
    bool isnotall = (toexcludefrom[regnum] != -1);
    
    if (isnotall)
        myphysicalregions->errorundefined({toexcludefrom[regnum]});
    myphysicalregions->errorundefined(toexclude[regnum]);
        
    // Make sure the regions are of same dimension:
    int physregdim;
    if (isnotall)
        physregdim = myphysicalregions->get(toexcludefrom[regnum])->getelementdimension();
    else
        physregdim = myelements->getdimension();
    
    for (int i = 0; i < toexclude[regnum].size(); i++)
    {
        int curdim = myphysicalregions->get(toexclude[regnum][i])->getelementdimension();
        if (curdim != physregdim)
        {
            logs log;
            log.msg() << "Error in 'regiondefiner' object: cannot exclude a " << curdim << "D region form a " << physregdim << "D region (dimensions must be equal)" << std::endl;
            log.error();
        }
    }
    
    physicalregion* newphysreg = myphysicalregions->get(excluded[regnum]);
    physicalregion* curphysreg;
    if (isnotall)
        curphysreg = myphysicalregions->get(toexcludefrom[regnum]);
    
    std::vector<std::vector<int>>* curelems;
    if (isnotall)
        curelems = curphysreg->getelementlist();

    // Loop on all element types:
    for (int i = 0; i <= 7; i++)
    {
        int numelems;
        if (isnotall)
            numelems = curelems->at(i).size();
        else
            numelems = myelements->countcells(i);
        
        if (numelems == 0)
            continue;

        int numelemsintype = myelements->count(i);
        std::vector<bool> inexcluded(numelemsintype, false);
        // First add all elements from which to exclude:
        for (int e = 0; e < numelems; e++)
        {
            if (isnotall)
                inexcluded[curelems->at(i)[e]] = true;
            else
                inexcluded[e] = true;
        }

        // Now remove the elements to exclude:
        for (int j = 0; j < toexclude[regnum].size(); j++)
        {
            physicalregion* curtoexclude = myphysicalregions->get(toexclude[regnum][j]);
            std::vector<int>* curelemtypetoexclude = &(curtoexclude->getelementlist()->at(i));

            for (int e = 0; e < curelemtypetoexclude->size(); e++)
                inexcluded[curelemtypetoexclude->at(e)] = false;
        }

        // All trues must be added:
        for (int j = 0; j < inexcluded.size(); j++)
        {
            if (inexcluded[j])
                newphysreg->addelement(i, j);
        }
    }
    newphysreg->removeduplicatedelements();
}

void regiondefiner::definelayerregion(int regnum)
{
    bool isnotall = (tolayer[regnum] != -1);
    
    if (isnotall)
        myphysicalregions->errorundefined({tolayer[regnum]});
    myphysicalregions->errorundefined({growthstart[regnum]});
        
    physicalregion* newphysreg = myphysicalregions->get(layered[regnum]);
    physicalregion* origphysreg;
    if (isnotall)
        origphysreg = myphysicalregions->get(tolayer[regnum]);
    physicalregion* growthphysreg = myphysicalregions->get(growthstart[regnum]);
    
    int nl = numlayers[regnum];
    
    // Tag the nodes that are in the growth region:
    std::vector<bool> isnodeingrowthregion(mynodes->count(), false);
    std::vector<std::vector<int>>* elemsingr = growthphysreg->getelementlist();
    // Loop on all element types:
    for (int i = 0; i <= 7; i++)
    {
        int numelems = elemsingr->at(i).size();
        
        if (numelems == 0)
            continue;
            
        element el(i);
        int nn = el.countnodes();
            
        for (int j = 0; j < numelems; j++)
        {
            int curelem = elemsingr->at(i)[j];
            for (int n = 0; n < nn; n++)
                isnodeingrowthregion[myelements->getsubelement(0, i, curelem, n)] = true;
        }
    }
    
    std::vector<std::vector<bool>> inlayer(8, std::vector<bool>(0)); // avoid duplicates
    for (int i = 0; i <= 7; i++)
        inlayer[i] = std::vector<bool>(myelements->count(i), false);

    std::vector<std::vector<int>>* curelems;
    if (isnotall)
        curelems = origphysreg->getelementlist();
    
    for (int l = 0; l < nl; l++)
    {
        // 'isnodeingr' is fixed:
        std::vector<bool> isnodeingr = isnodeingrowthregion;
    
        // Find all elements in the current layer (must be touching the growth region with at least one node):
        for (int i = 0; i <= 7; i++)
        {
            int numelems;
            if (isnotall)
                numelems = curelems->at(i).size();
            else
                numelems = myelements->countcells(i);
            
            if (numelems == 0)
                continue;

            element el(i);
            int nn = el.countnodes();

            for (int j = 0; j < numelems; j++)
            {
                int curelem;
                if (isnotall)
                    curelem = curelems->at(i)[j];
                else
                    curelem = j;
                    
                if (inlayer[i][curelem])
                    continue;
                    
                for (int n = 0; n < nn; n++)
                {
                    int curnode = myelements->getsubelement(0, i, curelem, n);
                    // If element is in layer:
                    if (isnodeingr[curnode])
                    {
                        inlayer[i][curelem] = true;
                        newphysreg->addelement(i, curelem);
                        break;
                    }
                }
                // Add nodes if element is in layer:
                if (inlayer[i][curelem])
                {
                    for (int n = 0; n < nn; n++)
                    {
                        int curnode = myelements->getsubelement(0, i, curelem, n);
                        isnodeingrowthregion[curnode] = true;
                    }
                }
            }
        }
    }
    newphysreg->removeduplicatedelements();
}

void regiondefiner::defineanynoderegion(int regnum)
{
    bool isnotall = (toanynode[regnum] != -1);
    
    if (isnotall)
        myphysicalregions->errorundefined({toanynode[regnum]});
        
    physicalregion* newphysreg = myphysicalregions->get(anynoded[regnum]);
    physicalregion* origphysreg;
    if (isnotall)
        origphysreg = myphysicalregions->get(toanynode[regnum]);
    
    std::vector<std::vector<int>>* elemsinor;
    if (isnotall)
        elemsinor = origphysreg->getelementlist();
    
    // Loop on all element types:
    for (int i = 0; i <= 7; i++)
    {
        int numelems;
        if (isnotall)
            numelems = elemsinor->at(i).size();
        else
            numelems = myelements->countcells(i);
        
        if (numelems > 0)
        {
            int firstelem;
            if (isnotall)
                firstelem = elemsinor->at(i)[0];
            else
                firstelem = 0;
        
            int selnode = myelements->getsubelement(0, i, firstelem, 0);
            newphysreg->addelement(0, selnode);
         
            break;   
        }
    }
    newphysreg->removeduplicatedelements();
}

regiondefiner::regiondefiner(nodes& inputnodes, elements& inputelems, physicalregions& inputphysregs)
{
    mynodes = &inputnodes;
    myelements = &inputelems;
    myphysicalregions = &inputphysregs;
}

void regiondefiner::regionskin(int newphysreg, int physregtoskin)
{    
    int cur = toskin.size();
    std::vector<int> prio = {0,cur};
    mypriority.push_back(prio);

    skins.push_back(newphysreg);
    toskin.push_back(physregtoskin);
}

void regiondefiner::regionbox(int newphysreg, int selecteddim, std::vector<double> boxlimit, int physregtobox)
{
    int cur = tobox.size();
    std::vector<int> prio = {1,cur};
    mypriority.push_back(prio);
    
    if (boxlimit.size() != 6)
    {
        logs log;
        log.msg() << "Error in 'regiondefiner' object: expected a vector of length 6 for the box limits {x1,x2,y1,y2,z1,z2}" << std::endl;
        log.error();
    }
    if (selecteddim > 3 || selecteddim < 0)
    {
        logs log;
        log.msg() << "Error in 'regiondefiner' object: dimension of the elements to select cannot be " << selecteddim << std::endl;
        log.error();
    }

    boxed.push_back(newphysreg);
    tobox.push_back(physregtobox);
    boxelemdims.push_back(selecteddim);

    boxlimits.push_back(boxlimit);
}

void regiondefiner::regionsphere(int newphysreg, int selecteddim, std::vector<double> centercoords, double radius, int physregtosphere)
{
    int cur = tosphere.size();
    std::vector<int> prio = {2,cur};
    mypriority.push_back(prio);
    
    if (centercoords.size() != 3)
    {
        logs log;
        log.msg() << "Error in 'regiondefiner' object: expected a vector of length 3 for the sphere center {xc,yc,zc}" << std::endl;
        log.error();
    }
    if (selecteddim > 3 || selecteddim < 0)
    {
        logs log;
        log.msg() << "Error in 'regiondefiner' object: dimension of the elements to select cannot be " << selecteddim << std::endl;
        log.error();
    }

    sphered.push_back(newphysreg);
    tosphere.push_back(physregtosphere);
    sphereelemdims.push_back(selecteddim);

    spherecenters.push_back(centercoords);
    sphereradii.push_back(radius);
}

void regiondefiner::regionexclusion(int newphysreg, int physregtoexcludefrom, std::vector<int> physregstoexclude)
{    
    int cur = toexcludefrom.size();
    std::vector<int> prio = {3,cur};
    mypriority.push_back(prio);
    
    excluded.push_back(newphysreg);
    toexcludefrom.push_back(physregtoexcludefrom);
    toexclude.push_back(physregstoexclude);
}

void regiondefiner::regionlayer(int newphysreg, int physregtoselectfrom, int physregtostartgrowth, int nl)
{
    int cur = tolayer.size();
    std::vector<int> prio = {4,cur};
    mypriority.push_back(prio);
    
    if (nl <= 0)
    {
        logs log;
        log.msg() << "Error in 'regiondefiner' object: expected at least one layer to select" << std::endl;
        log.error();
    }

    layered.push_back(newphysreg);
    tolayer.push_back(physregtoselectfrom);
    growthstart.push_back(physregtostartgrowth);
    numlayers.push_back(nl);
}

void regiondefiner::regionanynode(int newphysreg, int physregtoselectfrom)
{
    int cur = toanynode.size();
    std::vector<int> prio = {5,cur};
    mypriority.push_back(prio);
    
    anynoded.push_back(newphysreg);
    toanynode.push_back(physregtoselectfrom);
}

bool regiondefiner::isanyregiondefined(void)
{
    return (mypriority.size() > 0);
}

bool regiondefiner::isanycoordinatedependentregiondefined(void)
{
    for (int i = 0; i < mypriority.size(); i++)
    {
        if (mypriority[i][0] == 1 || mypriority[i][0] == 2)
            return true;
    }

    return false;
}


void regiondefiner::defineregions(void)
{
    if (isanycoordinatedependentregiondefined())
    {
        std::vector<double> nt = mynodes->getnoisethreshold();
        noisethreshold = gentools::sum(nt);
    }

    for (int i = 0; i < mypriority.size(); i++)
    {
        if (mypriority[i][0] == 0)
            defineskinregion(mypriority[i][1]);
        if (mypriority[i][0] == 1)
            defineboxregion(mypriority[i][1]);
        if (mypriority[i][0] == 2)
            definesphereregion(mypriority[i][1]);
        if (mypriority[i][0] == 3)
            defineexclusionregion(mypriority[i][1]);
        if (mypriority[i][0] == 4)
            definelayerregion(mypriority[i][1]);
        if (mypriority[i][0] == 5)
            defineanynoderegion(mypriority[i][1]);
    }
    
    noisethreshold = -1;
}

regiondefiner regiondefiner::copy(nodes* nds, elements* els, physicalregions* prs)
{
    regiondefiner out;
    
    out = *this;
    
    out.mynodes = nds;
    out.myelements = els;
    out.myphysicalregions = prs;
    
    return out;
}


