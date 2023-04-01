#include "iodata.h"


iodata::iodata(int interpolorder, int geointerpolorder, bool isitscalardata, std::vector<double> timevals)
{
    if (interpolorder <= 0 || geointerpolorder <= 0)
    {
        logs log;
        log.msg() << "Error in 'iodata' object: cannot have a negative or zero interpolation order" << std::endl;
        log.error();
    }

    myinterpolorder = interpolorder;
    mygeointerpolorder = geointerpolorder;
    isscalardata = isitscalardata;
    mytimevals = timevals;

     mycoords = std::vector<std::vector<std::vector<densemat>>>(3, std::vector<std::vector<densemat>>(8, std::vector<densemat>(0)));
     if (isscalardata)
        mydata = std::vector<std::vector<std::vector<densemat>>>(1, std::vector<std::vector<densemat>>(8, std::vector<densemat>(0)));
    else
        mydata = std::vector<std::vector<std::vector<densemat>>>(3, std::vector<std::vector<densemat>>(8, std::vector<densemat>(0)));
}

void iodata::combine(void)
{
    // Loop on every element type:
    for (int i = 0; i < 8; i++)
    {
        // Skip if there is a single block:
        if (mycoords[0][i].size() <= 1)
            continue;
    
        for (int s = 0; s < 3; s++)
            mycoords[s][i] = {densemat(mycoords[s][i])};
        for (int comp = 0; comp < mydata.size(); comp++)
            mydata[comp][i] = {densemat(mydata[comp][i])};
    }
}

bool iodata::isscalar(void) { return isscalardata; }
int iodata::getinterpolorder(void) { return myinterpolorder; };
int iodata::getgeointerpolorder(void) { return mygeointerpolorder; };

std::vector<double> iodata::gettimetags(void) { return mytimevals; };

bool iodata::ispopulated(int elemtypenum)
{
    return (mycoords[0][elemtypenum].size() > 0);
}

std::vector<int> iodata::getactiveelementtypes(void)
{
    std::vector<int> activeelementtypes = {};
    for (int i = 0; i < 8; i++)
    {
        if (ispopulated(i) == true)
            activeelementtypes.push_back(i);
    }
    return activeelementtypes;
}

void iodata::addcoordinates(int elemtypenum, densemat xcoords, densemat ycoords, densemat zcoords)
{
    element myelement(elemtypenum, mygeointerpolorder);
    int numgeonodes = myelement.countcurvednodes();
    
    numcoordtimestepsprovided  = xcoords.countcolumns()/numgeonodes;

    mycoords[0][elemtypenum].push_back(xcoords);
    mycoords[1][elemtypenum].push_back(ycoords);
    mycoords[2][elemtypenum].push_back(zcoords);
}

void iodata::adddata(int elemtypenum, std::vector<densemat> vals)
{
    element myelement(elemtypenum, myinterpolorder);
    int numdatanodes = myelement.countcurvednodes();
    
    numdatatimestepsprovided  = vals[0].countcolumns()/numdatanodes;
    
    // The non-provided components are set to 0:
    int vallen = vals.size();
    if (isscalardata == false && vallen < 3)
    {
        for (int i = vallen; i < 3; i++)
            vals.push_back(densemat(vals[0].countrows(), vals[0].countcolumns(), 0.0));
    }

    for (int comp = 0; comp < vals.size(); comp++)
        mydata[comp][elemtypenum].push_back(vals[comp]);
}

int iodata::countcoordnodes(int elemtypenum)
{
    combine();
    if (mycoords[0][elemtypenum].size() > 0)
        return mycoords[0][elemtypenum][0].count()/numcoordtimestepsprovided;
    else
        return 0;
}

int iodata::countcoordnodes(void)
{
    combine();
    int numnodes = 0;
    for (int i = 0; i < 8; i++)
        numnodes += countcoordnodes(i);
    return numnodes;
}

int iodata::countelements(int elemtypenum)
{
    combine();
    if (mycoords[0][elemtypenum].size() > 0)
        return mycoords[0][elemtypenum][0].countrows();
    else
        return 0;
}

int iodata::countelements(void)
{
    combine();
    int numelems = 0;
    for (int i = 0; i < 8; i++)
        numelems += countelements(i);
    return numelems;
}

std::vector<densemat> iodata::getcoordinates(int elemtypenum, int timestepindex)
{
    combine();
    
    if (timestepindex == -1 || numcoordtimestepsprovided == 1)
        return {mycoords[0][elemtypenum][0], mycoords[1][elemtypenum][0], mycoords[2][elemtypenum][0]};
    else
    {
        element myelement(elemtypenum, mygeointerpolorder);
        int numnod = myelement.countcurvednodes();
        int col1 = timestepindex*numnod, col2 = (timestepindex+1)*numnod-1;
        return {mycoords[0][elemtypenum][0].extractcols(col1,col2), mycoords[1][elemtypenum][0].extractcols(col1,col2), mycoords[2][elemtypenum][0].extractcols(col1,col2)};
    }
}

std::vector<densemat> iodata::getdata(int elemtypenum, int timestepindex)
{
    combine();
    
    if (timestepindex == -1 || numdatatimestepsprovided == 1)
    {
        if (isscalardata == true)
            return {mydata[0][elemtypenum][0]};
        else
            return {mydata[0][elemtypenum][0], mydata[1][elemtypenum][0], mydata[2][elemtypenum][0]};
    }
    else
    {
        element myelement(elemtypenum, myinterpolorder);
        int numnod = myelement.countcurvednodes();
        int col1 = timestepindex*numnod, col2 = (timestepindex+1)*numnod-1;
        
        if (isscalardata == true)
            return {mydata[0][elemtypenum][0].extractcols(col1,col2)};
        else
            return {mydata[0][elemtypenum][0].extractcols(col1,col2), mydata[1][elemtypenum][0].extractcols(col1,col2), mydata[2][elemtypenum][0].extractcols(col1,col2)};
    }
}

