#include "opon.h"


opon::opon(int physreg, expression* coordshift, std::shared_ptr<operation> arg)
{
    myphysreg = physreg; 
    myarg = arg;
    if (coordshift != NULL)
        mycoordshift = {*coordshift};
}

std::vector<std::vector<densematrix>> opon::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Otherwise the stored data will be used during the interpolation step (on wrong ki, eta, phi coordinates):
    universe::forbidreuse();
    

    // Calculate the x, y and z coordinates at which to interpolate:
    int problemdimension = universe::mymesh->getmeshdimension();
    field x("x"), y("y"), z("z");
    expression xdef, ydef, zdef;
    xdef = x; ydef = y; zdef = z;
    
    if (mycoordshift.size() > 0)
    {
        xdef = x+mathop::compx(mycoordshift[0]);
        if (meshdeform->countrows() > 1)
            ydef = y+mathop::compy(mycoordshift[0]);
        if (meshdeform->countrows() > 2)
            zdef = z+mathop::compz(mycoordshift[0]);
    }   
    expression xyz(3, 1, {xdef,ydef,zdef});
    
    densematrix xvalmat, yvalmat, zvalmat;
    
    xvalmat = xdef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL)[1][0];
    yvalmat = ydef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL)[1][0];
    zvalmat = zdef.getoperationinarray(0,0)->interpolate(elemselect, evaluationcoordinates, NULL)[1][0];
    
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
    std::vector<double> interpolated;
    std::vector<bool> isfound;
    
    if (meshdeform == NULL)
        expression(myarg).interpolate(myphysreg, xyzcoords, interpolated, isfound); 
    else
        expression(myarg).interpolate(myphysreg, *meshdeform, xyzcoords, interpolated, isfound); 
    
    for (int i = 0; i < isfound.size(); i++)
    {
        if (isfound[i] == false)
        {
            std::cout << "Error in 'opon' object: trying to interpolate at a point outside of physical region " << myphysreg << " (x,y,z) = (" << interpolated[3*i+0] << ", " << interpolated[3*i+1] << ", " << interpolated[3*i+2] << ")" << std::endl;
            abort();
        }
    }
    
    
    // Place the interpolated values in a densematrix:
    densematrix output(numelems, numgp);
    double* outputvals = output.getvalues();
    
    for (int i = 0; i < numpts; i++)
        outputvals[i] = interpolated[i];
    
    std::vector<std::vector<densematrix>> outvec = {{},{output}};
    
    
    universe::allowreuse();    
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), outvec);
        
    return outvec;
}

densematrix opon::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    std::cout << "Error in 'opon' object: cannot perform a multiharmonic interpolation" << std::endl;
    abort();
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

    return shared_from_this();
}

std::shared_ptr<operation> opon::copy(void)
{
    expression* exprptr = NULL;
    if (mycoordshift.size() > 0)
        exprptr = &(mycoordshift[0]);

    std::shared_ptr<opon> op(new opon(myphysreg, exprptr, myarg));
    *op = *this;
    op->reuse = false;
    return op;
}

void opon::print(void)
{
    std::cout << "on(" << myphysreg << ", ";
    if (mycoordshift.size() > 0)
        std::cout << "movedxyz, ";
    myarg->print();
    std::cout << ")";
}
