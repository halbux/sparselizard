#include "regiondefiner.h"


void regiondefiner::defineskinregion(int regnum)
{
    myphysicalregions->errorundefined({toskin[regnum]});

    physicalregion* newphysreg = myphysicalregions->get(skins[regnum]);
    physicalregion* curphysreg = myphysicalregions->get(toskin[regnum]);

    if (curphysreg->getelementdimension() == 0)
    {
        std::cout << "Error in 'regiondefiner' object: cannot get the skin of point elements" << std::endl;
        abort();
    }

    int physregdim = curphysreg->getelementdimension();
    // The skin has a one lower dimension:
    int skindim = physregdim-1;

    std::vector<std::vector<int>>* curelems = curphysreg->getelementlist();

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
            std::vector<int>* curelemtype = &(curelems->at(elemtype));

            int numelems = curelemtype->size();
            if (numelems == 0)
                continue;

            element myelement(elemtype);

            int numsubtypeelems = myelement.counttype(skinelemtype);
            if (numsubtypeelems == 0)
                continue;

            // Loop on all elements:
            for (int elem = 0; elem < numelems; elem++)
            {
                int curelem = curelemtype->at(elem);

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
    myphysicalregions->errorundefined({tobox[regnum]});

    std::vector<double>* nodecoords = mynodes->getcoordinates();

    std::vector<double> boxlimit = boxlimits[regnum];

    physicalregion* newphysreg = myphysicalregions->get(boxed[regnum]);
    physicalregion* curphysreg = myphysicalregions->get(tobox[regnum]);

    std::vector<std::vector<int>>* curelems = curphysreg->getelementlist();

    // Loop on all element types:
    for (int elemtype = 0; elemtype <= 7; elemtype++)
    {
        if (curelems->at(elemtype).size() == 0)
            continue;

        element myelement(elemtype);

        std::vector<int>* curelemtype = &(curelems->at(elemtype));

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
            for (int elem = 0; elem < curelemtype->size(); elem++)
            {
                int curelem = curelemtype->at(elem);

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
    myphysicalregions->errorundefined({tosphere[regnum]});

    std::vector<double>* nodecoords = mynodes->getcoordinates();

    std::vector<double> spherecenter = spherecenters[regnum];
    double sphereradius = sphereradii[regnum];

    physicalregion* newphysreg = myphysicalregions->get(sphered[regnum]);
    physicalregion* curphysreg = myphysicalregions->get(tosphere[regnum]);

    std::vector<std::vector<int>>* curelems = curphysreg->getelementlist();

    // Loop on all element types:
    for (int elemtype = 0; elemtype <= 7; elemtype++)
    {
        if (curelems->at(elemtype).size() == 0)
            continue;

        element myelement(elemtype);

        std::vector<int>* curelemtype = &(curelems->at(elemtype));

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
            for (int elem = 0; elem < curelemtype->size(); elem++)
            {
                int curelem = curelemtype->at(elem);

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

void regiondefiner::boxselection(int newphysreg, int selecteddim, std::vector<double> boxlimit, int physregtobox)
{
    int cur = tobox.size();
    std::vector<int> prio = {1,cur};
    mypriority.push_back(prio);
    
    if (boxlimit.size() != 6)
    {
        std::cout << "Error in 'regiondefiner' object: expected a vector of length 6 for the box limits {x1,x2,y1,y2,z1,z2}" << std::endl;
        abort();
    }
    if (selecteddim > 3 || selecteddim < 0)
    {
        std::cout << "Error in 'regiondefiner' object: dimension of the elements to select cannot be " << selecteddim << std::endl;
        abort();
    }

    boxed.push_back(newphysreg);
    tobox.push_back(physregtobox);
    boxelemdims.push_back(selecteddim);

    // Make the box limit slightly larger to remove the roundoff noise issues:
    boxlimit[0] -= boxlimit[0]*roundoffnoise; boxlimit[2] -= boxlimit[2]*roundoffnoise; boxlimit[4] -= boxlimit[4]*roundoffnoise;
    boxlimit[1] += boxlimit[1]*roundoffnoise; boxlimit[3] += boxlimit[3]*roundoffnoise; boxlimit[5] += boxlimit[5]*roundoffnoise;

    boxlimits.push_back(boxlimit);
}

void regiondefiner::sphereselection(int newphysreg, int selecteddim, std::vector<double> centercoords, double radius, int physregtosphere)
{
    int cur = tosphere.size();
    std::vector<int> prio = {2,cur};
    mypriority.push_back(prio);
    
    if (centercoords.size() != 3)
    {
        std::cout << "Error in 'regiondefiner' object: expected a vector of length 3 for the sphere center {xc,yc,zc}" << std::endl;
        abort();
    }
    if (selecteddim > 3 || selecteddim < 0)
    {
        std::cout << "Error in 'regiondefiner' object: dimension of the elements to select cannot be " << selecteddim << std::endl;
        abort();
    }

    sphered.push_back(newphysreg);
    tosphere.push_back(physregtosphere);
    sphereelemdims.push_back(selecteddim);

    // Make the sphere radius slightly larger to remove the roundoff noise issues:
    radius += radius*roundoffnoise;

    spherecenters.push_back(centercoords);
    sphereradii.push_back(radius);
}


void regiondefiner::defineregions(void)
{
    for (int i = 0; i < mypriority.size(); i++)
    {
        if (mypriority[i][0] == 0)
            defineskinregion(mypriority[i][1]);
        if (mypriority[i][0] == 1)
            defineboxregion(mypriority[i][1]);
        if (mypriority[i][0] == 2)
            definesphereregion(mypriority[i][1]);
    }
}


