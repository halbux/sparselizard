#include "opon.h"


opon::opon(int physreg, expression* coordshift, std::shared_ptr<operation> arg, bool errorifnotfound)
{
    myphysreg = physreg; 
    myerrorifnotfound = errorifnotfound;
    myarg = arg;
    if (coordshift != NULL)
        mycoordshift = {*coordshift};
}

std::vector<std::vector<densemat>> opon::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }

    // Otherwise the stored data will be used during the interpolation step (on wrong ki, eta, phi coordinates):
    bool wasreuseallowed = universe::isreuseallowed;
    universe::isreuseallowed = false;
    

    // Calculate the x, y and z coordinates at which to interpolate:
    field x("x"), y("y"), z("z");
    expression xdef, ydef, zdef;
    xdef = x; ydef = y; zdef = z;
    
    if (mycoordshift.size() > 0)
    {
        xdef = x+sl::compx(mycoordshift[0]);
        if (mycoordshift[0].countrows() > 1)
            ydef = y+sl::compy(mycoordshift[0]);
        if (mycoordshift[0].countrows() > 2)
            zdef = z+sl::compz(mycoordshift[0]);
    }   
    expression xyz(3, 1, {xdef,ydef,zdef});
    
    std::vector<std::vector<densemat>> xmatvec, ymatvec, zmatvec;
    
    xmatvec = xdef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL);
    ymatvec = ydef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL);
    zmatvec = zdef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL);
    
    // Make sure the coordinate shift expression was constant in time:
    if (xmatvec.size() > 2 || ymatvec.size() > 2 || zmatvec.size() > 2)
    {
        logs log;
        log.msg() << "Error in 'opon' object: coordinate shift expression cannot be (multi) harmonic" << std::endl;
        log.error();
    }
    
    densemat xvalmat, yvalmat, zvalmat;
    xvalmat = xmatvec[1][0];
    yvalmat = ymatvec[1][0];
    zvalmat = zmatvec[1][0];
    
    double* xvals = xvalmat.getvalues();
    double* yvals = yvalmat.getvalues();
    double* zvals = zvalmat.getvalues();
    
    
    int numelems = xvalmat.countrows(), numgp = xvalmat.countcolumns();
    int numpts = numelems*numgp;
    
    std::vector<double> xyzcoords(3*numpts);
    for (int i = 0; i < numpts; i++)
    {
        xyzcoords[3*i+0] = xvals[i];
        xyzcoords[3*i+1] = yvals[i];
        xyzcoords[3*i+2] = zvals[i];
    }
    
    
    // Interpolate the expression at these coordinates:
    std::vector<std::vector<double>> interpolated;
    std::vector<bool> isfound;
    
    expression(myarg).interpolate(myphysreg, meshdeform, xyzcoords, interpolated, isfound, -1); 
    
    if (myerrorifnotfound)
    {
        for (int i = 0; i < isfound.size(); i++)
        {
            if (isfound[i] == false)
            {
                logs log;
                log.msg() << "Error in 'opon' object: trying to interpolate at a point outside of physical region " << myphysreg << " or interpolation algorithm failed to converge" << std::endl;
                log.msg() << "Error was at (x,y,z) = (" << xyzcoords[3*i+0] << ", " << xyzcoords[3*i+1] << ", " << xyzcoords[3*i+2] << ")" << std::endl;
                log.error();
            }
        }
    }
    
    
    // Place the interpolated values in a std::vector<std::vector<densemat>>:
    std::vector<std::vector<densemat>> outvec(interpolated.size(), std::vector<densemat>(0));
    for (int h = 0; h < interpolated.size(); h++)
    {
        if (interpolated[h].size() > 0)
        {
            densemat singleharm(numelems, numgp);
            double* singleharmvals = singleharm.getvalues();

            for (int i = 0; i < numpts; i++)
                singleharmvals[i] = interpolated[h][i];
            outvec[h] = {singleharm};
        }
    }
    
    
    universe::isreuseallowed = wasreuseallowed;   
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), outvec);
        
    return outvec;
}

densemat opon::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvaluefft(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputedfft(precomputedindex); }
    }

    // Otherwise the stored data will be used during the interpolation step (on wrong ki, eta, phi coordinates):
    bool wasreuseallowed = universe::isreuseallowed;
    universe::isreuseallowed = false;
    

    // Calculate the x, y and z coordinates at which to interpolate:
    field x("x"), y("y"), z("z");
    expression xdef, ydef, zdef;
    xdef = x; ydef = y; zdef = z;
    
    if (mycoordshift.size() > 0)
    {
        xdef = x+sl::compx(mycoordshift[0]);
        if (mycoordshift[0].countrows() > 1)
            ydef = y+sl::compy(mycoordshift[0]);
        if (mycoordshift[0].countrows() > 2)
            zdef = z+sl::compz(mycoordshift[0]);
    }   
    expression xyz(3, 1, {xdef,ydef,zdef});
    
    std::vector<std::vector<densemat>> xmatvec, ymatvec, zmatvec;
    
    xmatvec = xdef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL);
    ymatvec = ydef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL);
    zmatvec = zdef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL);
    
    // Make sure the coordinate shift expression was constant in time:
    if (xmatvec.size() > 2 || ymatvec.size() > 2 || zmatvec.size() > 2)
    {
        logs log;
        log.msg() << "Error in 'opon' object: coordinate shift expression cannot be (multi) harmonic" << std::endl;
        log.error();
    }
    
    densemat xvalmat, yvalmat, zvalmat;
    xvalmat = xmatvec[1][0];
    yvalmat = ymatvec[1][0];
    zvalmat = zmatvec[1][0];
    
    double* xvals = xvalmat.getvalues();
    double* yvals = yvalmat.getvalues();
    double* zvals = zvalmat.getvalues();
    
    
    int numelems = xvalmat.countrows(), numgp = xvalmat.countcolumns();
    int numpts = numelems*numgp;
    
    std::vector<double> xyzcoords(3*numpts);
    for (int i = 0; i < numpts; i++)
    {
        xyzcoords[3*i+0] = xvals[i];
        xyzcoords[3*i+1] = yvals[i];
        xyzcoords[3*i+2] = zvals[i];
    }
    
    
    // Interpolate the expression at these coordinates:
    std::vector<std::vector<double>> interpolated;
    std::vector<bool> isfound;
    
    expression(myarg).interpolate(myphysreg, meshdeform, xyzcoords, interpolated, isfound, numtimeevals);

    if (myerrorifnotfound)
    {
        for (int i = 0; i < isfound.size(); i++)
        {
            if (isfound[i] == false)
            {
                logs log;
                log.msg() << "Error in 'opon' object: trying to interpolate at a point outside of physical region " << myphysreg << " or interpolation algorithm failed to converge" << std::endl;
                log.msg() << "Error was at (x,y,z) = (" << xyzcoords[3*i+0] << ", " << xyzcoords[3*i+1] << ", " << xyzcoords[3*i+2] << ")" << std::endl;
                log.error();
            }
        }
    }
    
    // Place the interpolated values in a densemat:
    densemat outmat(numtimeevals, numpts);
    double* outmatvals = outmat.getvalues();
    
    for (int i = 0; i < numpts*numtimeevals; i++)
        outmatvals[i] = interpolated[0][i];
    
    
    universe::isreuseallowed = wasreuseallowed;
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputedfft(shared_from_this(), outmat);
        
    return outmat;
}

std::vector<std::shared_ptr<operation>> opon::getarguments(void)
{
    int numargs = 1;
    if (mycoordshift.size() > 0)
        numargs += mycoordshift[0].countrows();
    std::vector<std::shared_ptr<operation>> allargs(numargs);
    allargs[0] = myarg;
    if (mycoordshift.size() > 0)
    {
        for (int i = 0; i < mycoordshift[0].countrows(); i++)
            allargs[1+i] = mycoordshift[0].getoperationinarray(i,0);
    }
    
    return allargs;
}

std::shared_ptr<operation> opon::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);

    if (myarg->isconstant())
        return std::shared_ptr<operation>(new opconstant(myarg->getvalue()));
    else
        return shared_from_this();
}

std::shared_ptr<operation> opon::copy(void)
{
    expression* exprptr = NULL;
    if (mycoordshift.size() > 0)
        exprptr = &(mycoordshift[0]);

    std::shared_ptr<opon> op(new opon(myphysreg, exprptr, myarg, myerrorifnotfound));
    *op = *this;
    op->reuse = false;
    return op;
}

void opon::print(void)
{
    std::cout << "on(" << myphysreg << ", ";
    myarg->print();
    std::cout << ")";
}
