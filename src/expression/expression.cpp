#include "expression.h"
#include "oncontext.h"


expression::expression(field input)
{
    mynumrows = input.countcomponents();
    mynumcols = 1;

    myoperations.resize(mynumrows);
    for (int row = 0; row < mynumrows; row++)
    {
        // For fields with multiple form function components:
        if (input.getpointer()->countformfunctioncomponents() > 1)
        {
            std::shared_ptr<opfield> op(new opfield(input.getpointer()));
            op->selectformfunctioncomponent(row);
            myoperations[row] = op;
        }
        else
        {
            std::shared_ptr<opfield> op(new opfield(input.comp(row).getpointer()));
            // For printing purposes:
            if (mynumrows > 1)
                op->setfieldcomponent(row);
            myoperations[row] = op;
        }
    }

    // Transform from the reference element to the physical one a hcurl field:
    if (input.getpointer()->gettypename() == "hcurl")
    {
        inrefcoord = {std::make_pair("hcurl", *this)};
        myoperations = ( invjac()*(*this) ).myoperations;
    }
}

expression::expression(double input) { mynumrows = 1; mynumcols = 1; myoperations = {std::shared_ptr<opconstant>(new opconstant(input))}; }

expression::expression(parameter& input)
{
    mynumrows = input.countrows();
    mynumcols = input.countcolumns();

    myoperations.resize(mynumrows*mynumcols);
    for (int row = 0; row < mynumrows; row++)
    {
        for (int col = 0; col < mynumcols; col++)
            myoperations[row*mynumcols+col] = std::shared_ptr<opparameter>(new opparameter(&input, row, col));
    }
}

expression::expression(int numrows, int numcols, std::vector<expression> input)
{
    mynumrows = numrows;
    mynumcols = numcols;

    // In case the user provides only the diagonal part of a square diagonal matrix:
    if (mynumrows == mynumcols && input.size() == mynumrows)
    {
        std::vector<expression> fullinput(mynumrows*mynumcols);
        for (int i = 0; i < mynumrows; i++)
        {
            for (int j = 0; j < mynumcols; j++)
            {
                if (i == j)
                    fullinput[i*mynumcols+j] = input[i];
                else
                    fullinput[i*mynumcols+j] = expression(0);
            }
        }
        input = fullinput;
    }

    // In case the user provides only the lower triangular part of a square symmetric matrix:
    if (mynumrows == mynumcols && input.size() == mynumrows*(mynumcols+1)/2)
    {
        std::vector<expression> fullinput(mynumrows*mynumcols);
        int index = 0;
        for (int i = 0; i < mynumrows; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (i == j)
                    fullinput[i*mynumcols+j] = input[index];
                else
                {
                    fullinput[i*mynumcols+j] = input[index];
                    fullinput[j*mynumcols+i] = input[index];
                }

                index++;
            }
        }
        input = fullinput;
    }

    if (mynumrows*mynumcols != input.size())
    {
        std::cout << "Error in 'expression' object: a vector of length " << mynumrows*mynumcols << " is required" << std::endl;
        abort();
    }

    myoperations.resize(mynumrows*mynumcols);
    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (input[i].isscalar())
            myoperations[i] = input[i].myoperations[0];
        else
        {
            std::cout << "Error in 'expression' object: expressions provided in vector must all be scalar" << std::endl;
            abort();
        }
    }
}

expression::expression(const std::vector<std::vector<expression>> input)
{
    if (input.size() == 0)
        return;

    // Single column special case:
    if (input.size() == 1)
    {
        if (input[0].size() == 0)
            return;

        // Get the final expression dimension:
        int numcols = input[0][0].mynumcols;
        int numrows = 0;
        for (int i = 0; i < input[0].size(); i++)
        {
            if (input[0][i].mynumcols != numcols)
            {
                std::cout << "Error in 'expression' object: expression dimension mismatch in concatenation" << std::endl;
                abort();
            }
            numrows += input[0][i].mynumrows;
        }

        mynumrows = numrows;
        mynumcols = numcols;
        myoperations.resize(numrows*numcols);

        int index = 0;
        for (int i = 0; i < input[0].size(); i++)
        {
            for (int j = 0; j < input[0][i].mynumrows*input[0][i].mynumcols; j++)
            {
                myoperations[index] = input[0][i].myoperations[j];
                index++;
            }
        }
        return;
    }

    // Combine all rows in each column:
    std::vector<std::vector<expression>> exprs(1,std::vector<expression>(input.size()));
    for (int i = 0; i < input.size(); i++)
    {
        expression rowsconcatenated({input[i]});
        exprs[0][i] = rowsconcatenated.transpose();
    }

    expression finalconcatenated(exprs);
    finalconcatenated = finalconcatenated.transpose();

    mynumrows = finalconcatenated.mynumrows;
    mynumcols = finalconcatenated.mynumcols;
    myoperations = finalconcatenated.myoperations;
}

expression::expression(expression condexpr, expression exprtrue, expression exprfalse)
{
    // Make sure the conditional expression is a scalar:
    if (condexpr.countrows() != 1 || condexpr.countcolumns() != 1)
    {
        std::cout << "Error in 'expression' object: expected a scalar condition expression (first argument)" << std::endl;
        abort();
    }
    // Make sure the two other expressions have the same size:
    if (exprtrue.countrows() != exprfalse.countrows() || exprtrue.countcolumns() != exprfalse.countcolumns())
    {
        std::cout << "Error in 'expression' object: expected same sized expressions in last two arguments" << std::endl;
        abort();
    }

    if (condexpr.myoperations[0]->isdofincluded() || condexpr.myoperations[0]->istfincluded())
    {
        std::cout << "Error in 'expression' object: conditional expression arguments cannot include a dof() or tf()" << std::endl;
        abort();
    }

    // Break into scalars:
    mynumrows = exprtrue.countrows(); mynumcols = exprtrue.countcolumns();

    myoperations.resize(exprtrue.countrows()*exprtrue.countcolumns());
    for (int i = 0; i < exprtrue.countrows()*exprtrue.countcolumns(); i++)
    {
        if (exprtrue.myoperations[i]->isdofincluded() || exprtrue.myoperations[i]->istfincluded() || exprfalse.myoperations[i]->isdofincluded() || exprfalse.myoperations[i]->istfincluded())
        {
            std::cout << "Error in 'expression' object: conditional expression arguments cannot include a dof() or tf()" << std::endl;
            abort();
        }
        myoperations[i] = std::shared_ptr<opcondition>(new opcondition(condexpr.myoperations[0], exprtrue.myoperations[i], exprfalse.myoperations[i]));
    }
}

expression::expression(spline spl, expression arg)
{
    if (arg.isscalar() == false)
    {
        std::cout << "Error in 'expression' object: expected a scalar expression as argument for the spline interpolation" << std::endl;
        abort();
    }
    if (arg.myoperations[0]->isdofincluded() || arg.myoperations[0]->istfincluded())
    {
        std::cout << "Error in 'expression' object: spline argument cannot include a dof() or tf()" << std::endl;
        abort();
    }
    mynumrows = 1; mynumcols = 1;
    myoperations = {std::shared_ptr<opspline>(new opspline(spl,arg.myoperations[0]))};
}

expression::expression(std::vector<double> pos, std::vector<expression> exprs, expression tocompare)
{
    int numintervals = pos.size()+1;
    if (numintervals != exprs.size())
    {
        std::cout << "Error in 'expression' object: expected " << numintervals << " expressions to define the " << numintervals << " intervals from -infinity to +infinity" << std::endl;
        abort();
    }
    if (tocompare.isscalar() == false)
    {
        std::cout << "Error in 'expression' object: expected a scalar expression as interval variable" << std::endl;
        abort();
    }
    if (tocompare.myoperations[0]->isdofincluded() || tocompare.myoperations[0]->istfincluded())
    {
        std::cout << "Error in 'expression' object: interval variable cannot include a dof() or tf()" << std::endl;
        abort();
    }
    // Make sure the positions are sorted ascendingly:
    for (int i = 1; i < pos.size(); i++)
    {
        if (pos[i] <= pos[i-1])
        {
            std::cout << "Error in 'expression' object: positions for intervals must be sorted ascendingly" << std::endl;
            abort(); 
        }
    }
    
    expression expr = exprs[0];
    
    tocompare.reuseit();
    if (numintervals > 1)
    {
        std::vector<double> posoneless = pos; posoneless.pop_back();
        std::vector<expression> exprsoneless = exprs; exprsoneless.pop_back();
        
        expression oneless(posoneless, exprsoneless, tocompare);
        
        expr = mathop::ifpositive(tocompare-pos[numintervals-2], exprs[numintervals-1], oneless);
    }
    
    mynumrows = expr.mynumrows; mynumcols = expr.mynumcols;
    myoperations = expr.myoperations;
}

expression::expression(int m, int n, std::vector<densematrix> customfct(std::vector<densematrix>), std::vector<expression> exprs)
{
    mynumrows = m; mynumcols = n;

    if (m <= 0 || n <= 0)
    {
        std::cout << "Error in 'expression' object: custom expression size cannot be zero or negative" << std::endl;
        abort();
    }
    
    // Get all argument operations:
    std::vector<std::shared_ptr<operation>> argops = {};
    for (int i = 0; i < exprs.size(); i++)
    {
        for (int j = 0; j < exprs[i].myoperations.size(); j++)
            argops.push_back(exprs[i].myoperations[j]);
    }
    
    myoperations.resize(m*n);
    std::vector<std::weak_ptr<opcustom>> weakops(m*n);
    for (int i = 0; i < m*n; i++)
    {
        std::shared_ptr<opcustom> op(new opcustom(i, customfct, argops));
        myoperations[i] = op;
        weakops[i] = op;
    }
    for (int i = 0; i < m*n; i++)
        (weakops[i].lock())->setfamily(weakops);
}

expression::expression(std::shared_ptr<operation> input)
{
    mynumrows = 1; mynumcols = 1;
    myoperations = {input};
}

expression expression::getrow(int rownum)
{
    if (rownum < 0 || rownum >= mynumrows)
    {
        std::cout << "Error in 'expression' object: cannot get row " << rownum << " in a " << mynumrows << "x" << mynumcols << " expression" << std::endl;
        abort();
    }

    expression output;
    output.mynumrows = 1;
    output.mynumcols = mynumcols;
    output.myoperations.resize(mynumcols);

    for (int i = 0; i < mynumcols; i++)
        output.myoperations[i] = myoperations[rownum*mynumcols+i];

    return output;
}

expression expression::getcolumn(int colnum)
{
    if (colnum < 0 || colnum >= mynumcols)
    {
        std::cout << "Error in 'expression' object: cannot get column " << colnum << " in a " << mynumrows << "x" << mynumcols << " expression" << std::endl;
        abort();
    }

    expression output;
    output.mynumrows = mynumrows;
    output.mynumcols = 1;
    output.myoperations.resize(mynumrows);

    for (int i = 0; i < mynumrows; i++)
        output.myoperations[i] = myoperations[i*mynumcols+colnum];

    return output;
}

void expression::reorderrows(std::vector<int> neworder)
{
    if (mynumrows != neworder.size())
    {
        std::cout << "Error in 'expression' object: cannot reorder rows with the vector provided (incorrect vector size)" << std::endl;
        abort();
    }

    int minval = *min_element(neworder.begin(), neworder.end());
    int maxval = *max_element(neworder.begin(), neworder.end());

    if (minval < 0 || maxval >= neworder.size())
    {
        std::cout << "Error in 'expression' object: cannot reorder rows with the vector provided (out of range integers)" << std::endl;
        abort();
    }

    std::vector<std::shared_ptr<operation>> ops = myoperations;

    for (int i = 0; i < mynumrows; i++)
    {
        for (int j = 0; j < mynumcols; j++)
            myoperations[i*mynumcols+j] = ops[neworder[i]*mynumcols+j];
    }
}

void expression::reordercolumns(std::vector<int> neworder)
{
    if (mynumcols != neworder.size())
    {
        std::cout << "Error in 'expression' object: cannot reorder columns with the vector provided (incorrect vector size)" << std::endl;
        abort();
    }

    int minval = *min_element(neworder.begin(), neworder.end());
    int maxval = *max_element(neworder.begin(), neworder.end());

    if (minval < 0 || maxval >= neworder.size())
    {
        std::cout << "Error in 'expression' object: cannot reorder columns with the vector provided (out of range integers)" << std::endl;
        abort();
    }

    std::vector<std::shared_ptr<operation>> ops = myoperations;

    for (int i = 0; i < mynumrows; i++)
    {
        for (int j = 0; j < mynumcols; j++)
            myoperations[i*mynumcols+j] = ops[i*mynumcols+neworder[j]];
    }
}


std::vector<double> expression::max(int physreg, int refinement, std::vector<double> xyzrange)
{
    return max(physreg, NULL, refinement, xyzrange);
}

std::vector<double> expression::max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{
    return max(physreg, &meshdeform, refinement, xyzrange);
}

std::vector<double> expression::min(int physreg, int refinement, std::vector<double> xyzrange)
{
    // The actual min value is minus the max value found:
    std::vector<double> output = (-this->getcopy()).max(physreg, NULL, refinement, xyzrange);
    if (output.size() > 0)
        output[0] = -output[0];
    return output;
}

std::vector<double> expression::min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{
    // The actual min value is minus the max value found:
    std::vector<double> output =  (-this->getcopy()).max(physreg, &meshdeform, refinement, xyzrange);
    if (output.size() > 0)
        output[0] = -output[0];
    return output;
}

std::vector<double> expression::max(int physreg, expression* meshdeform, int refinement, std::vector<double> xyzrange)
{
    field x("x"), y("y"), z("z");

    // Minimum refinement order is 1!
    if (refinement < 1)
        refinement = 1;

    // Make sure this expression is scalar and the mesh
    // deformation expression has the right size.
    if (not(isscalar()))
    {
        std::cout << "Error in 'expression' object: cannot get the max/min of a nonscalar expression" << std::endl;
        abort();
    }
    int problemdimension = universe::mymesh->getmeshdimension();
    if (meshdeform != NULL && (meshdeform->countcolumns() != 1 || meshdeform->countrows() < problemdimension))
    {
        std::cout << "Error in 'expression' object: mesh deformation expression has size " << meshdeform->countrows() << "x" << meshdeform->countcolumns() << " (expected " << problemdimension << "x1)" << std::endl;
        abort();
    }

    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();

    // Multiharmonic expressions are not allowed.
    if (not(isharmonicone(selecteddisjregs)))
    {
        std::cout << "Error in 'expression' object: cannot get the max/min of a multiharmonic expression (only constant harmonic 1)" << std::endl;
        abort();
    }
    if (meshdeform != NULL && not(meshdeform->isharmonicone(selecteddisjregs)))
    {
        std::cout << "Error in 'expression' object: the mesh deformation expression cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        abort();
    }
    
    universe::allowestimatorupdate(true);

    // This will be the output:
    std::vector<double> maxval = {};

    // Send the disjoint regions with same element type numbers together:
    disjointregionselector mydisjregselector(selecteddisjregs, {});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);

        // Get the node coordinates in the refined element:
        int elementtypenumber = (universe::mymesh->getdisjointregions())->getelementtypenumber(mydisjregs[0]);
        element myelement(elementtypenumber, refinement);
        std::vector<double> evaluationpoints = myelement.listnodecoordinates();

        // Loop on all total orientations (if required):
        bool isorientationdependent = isvalueorientationdependent(mydisjregs) || (meshdeform != NULL && meshdeform->isvalueorientationdependent(mydisjregs));
        elementselector myselector(mydisjregs, isorientationdependent);
        do
        {
            universe::allowreuse();

            // Compute the expression at the evaluation points:
            densematrix compxval = myoperations[0]->interpolate(myselector, evaluationpoints, meshdeform)[1][0];
            // Get the coordinates corresponding to the interpolated values:
            densematrix xval = expression(x).myoperations[0]->interpolate(myselector, evaluationpoints, meshdeform)[1][0];
            densematrix yval = expression(y).myoperations[0]->interpolate(myselector, evaluationpoints, meshdeform)[1][0];
            densematrix zval = expression(z).myoperations[0]->interpolate(myselector, evaluationpoints, meshdeform)[1][0];

            universe::forbidreuse();

            double* valuesptr = compxval.getvalues();
            double* xvalptr = xval.getvalues();
            double* yvalptr = yval.getvalues();
            double* zvalptr = zval.getvalues();

            // Loop on all data points.
            for (int d = 0; d < compxval.count(); d++)
            {
                bool isinboundedregion = (xyzrange.size() == 0) || (xyzrange[0] < xvalptr[d] && xyzrange[1] > xvalptr[d] && xyzrange[2] < yvalptr[d] && xyzrange[3] > yvalptr[d] && xyzrange[4] < zvalptr[d] && xyzrange[5] > zvalptr[d]);

                if ( isinboundedregion && (maxval.size() == 0 || maxval[0] < valuesptr[d]) )
                        maxval = {valuesptr[d], xvalptr[d], yvalptr[d], zvalptr[d]};
            }
        }
        while (myselector.next());
    }
    
    universe::allowestimatorupdate(false);

    return maxval;
}

void expression::interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{
    interpolate(physreg, NULL, xyzcoord, interpolated, isfound);
}

void expression::interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{
    interpolate(physreg, &meshdeform, xyzcoord, interpolated, isfound);
}

std::vector<double> expression::interpolate(int physreg, const std::vector<double> xyzcoord)
{
    std::vector<double> xyz = xyzcoord;

    if (xyz.size() != 3)
    {
        std::cout << "Error in 'expression' object: interpolate expected a coordinate vector of length 3" << std::endl;
        abort();
    }

    std::vector<double> interpolated = {};
    std::vector<bool> isfound = {};

    interpolate(physreg, NULL, xyz, interpolated, isfound);

    if (isfound[0])
        return interpolated;
    else
        return {};
}

std::vector<double> expression::interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord)
{
    std::vector<double> xyz = xyzcoord;

    if (xyz.size() != 3)
    {
        std::cout << "Error in 'expression' object: interpolate expected a coordinate vector of length 3" << std::endl;
        abort();
    }

    std::vector<double> interpolated = {};
    std::vector<bool> isfound = {};

    interpolate(physreg, &meshdeform, xyz, interpolated, isfound);

    if (isfound[0])
        return interpolated;
    else
        return {};
}

void expression::interpolate(int physreg, expression* meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{
    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> disjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();

    // Multiharmonic expressions are not allowed.
    if (not(isharmonicone(disjregs)))
    {
        std::cout << "Error in 'expression' object: cannot interpolate a multiharmonic expression (only constant harmonic 1)" << std::endl;
        abort();
    }
    if (xyzcoord.size()%3 != 0)
    {
        std::cout << "Error in 'expression' object: the interpolation coordinates vector should have a length that is a multiple of 3" << std::endl;
        abort();
    }

    int numcoords = xyzcoord.size()/3;
    int exprlen = countrows()*countcolumns();
    interpolated.resize(numcoords*exprlen);

    // Interpolate every expression entry:
    for (int i = 0; i < exprlen; i++)
    {
        std::vector<std::vector<double>> interpolatedscalar;
        expression(myoperations[i]).interpolate(physreg, meshdeform, xyzcoord, interpolatedscalar, isfound, -1);

        for (int j = 0; j < numcoords; j++)
            interpolated[j*exprlen+i] = interpolatedscalar[1][j];
    }
}


void expression::interpolate(int physreg, expression* meshdeform, std::vector<double>& xyzcoord, std::vector<std::vector<double>>& interpolated, std::vector<bool>& isfound, int numtimeevals)
{
    // Make sure the mesh deformation expression has the right size.
    int problemdimension = universe::mymesh->getmeshdimension();
    if (meshdeform != NULL && (meshdeform->countcolumns() != 1 || meshdeform->countrows() < problemdimension))
    {
        std::cout << "Error in 'expression' object: mesh deformation expression has size " << meshdeform->countrows() << "x" << meshdeform->countcolumns() << " (expected " << problemdimension << "x1)" << std::endl;
        abort();
    }

    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> disjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();
    if (meshdeform != NULL && not(meshdeform->isharmonicone(disjregs)))
    {
        std::cout << "Error in 'expression' object: the mesh deformation expression cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        abort();
    }


    // Preallocate the input containers:
    int numcoords = xyzcoord.size()/3;
    isfound = std::vector<bool>(numcoords,false);

    interpolated = {};
    if (numtimeevals != -1)
        interpolated = {std::vector<double>(numcoords*numtimeevals,0.0)};
    
    universe::allowestimatorupdate(true);
    
    referencecoordinategroup rcg(xyzcoord);
    
    // Send all disjoint regions with same element type number together:
    disjointregionselector mydisjregselector(disjregs, {});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> curdisjregs = mydisjregselector.getgroup(i);
        
        rcg.evalat(curdisjregs);
    
        // Simplify all coeffs for faster computation later on.
        // Also check if orientation matters.
        myoperations[0] = myoperations[0]->simplify(curdisjregs);
        bool isorientationdependent = (myoperations[0]->isvalueorientationdependent(curdisjregs) || (meshdeform != NULL && meshdeform->isvalueorientationdependent(curdisjregs)));
        
        while (rcg.next())
        {
            std::vector<double> kietaphi = rcg.getreferencecoordinates();
            std::vector<int> coordindexes = rcg.getcoordinatenumber();
            std::vector<int> elemens = rcg.getelements();
            int numrefcoords = kietaphi.size()/3;
            
            for (int c = 0; c < coordindexes.size(); c++)
                isfound[coordindexes[c]] = true;
        
            // Loop on all total orientations (if required):
            elementselector myselector(curdisjregs, elemens, isorientationdependent);
            do 
            {
                std::vector<int> origindexes = myselector.getoriginalindexes();
            
                if (numtimeevals == -1)
                {
                    // Clean storage before allowing reuse:
                    universe::forbidreuse();
                    universe::allowreuse();
                    std::vector<std::vector<densematrix>> interp = myoperations[0]->interpolate(myselector, kietaphi, meshdeform);
                    universe::forbidreuse();

                    if (interpolated.size() < interp.size())
                        interpolated.resize(interp.size());
                    
                    for (int h = 0; h < interp.size(); h++)
                    {
                        if (interp[h].size() > 0)
                        {
                            if (interpolated[h].size() == 0)
                                interpolated[h] = std::vector<double>(numcoords,0.0);
                                
                            for (int e = 0; e < origindexes.size(); e++)
                            {
                                for (int c = 0; c < numrefcoords; c++)
                                    interpolated[h][coordindexes[origindexes[e]*numrefcoords+c]] = interp[h][0].getvalues()[e*numrefcoords+c];
                            }
                        }
                    }
                }
                else
                {
                    // Clean storage before allowing reuse:
                    universe::forbidreuse();
                    universe::allowreuse();
                    densematrix interp = myoperations[0]->multiharmonicinterpolate(numtimeevals, myselector, kietaphi, meshdeform);
                    universe::forbidreuse();

                    for (int tim = 0; tim < numtimeevals; tim++)
                    {
                        for (int e = 0; e < origindexes.size(); e++)
                        {
                            for (int c = 0; c < numrefcoords; c++)
                                interpolated[0][numcoords*tim+coordindexes[origindexes[e]*numrefcoords+c]] = interp.getvalues()[origindexes.size()*numrefcoords*tim+e*numrefcoords+c];
                        }
                    }
                }
            }
            while (myselector.next());   
        }
    }
    
    universe::allowestimatorupdate(false);
    
    // Provide a non-empty interpolation vector even in case nothing was found:
    if (interpolated.size() == 0)
        interpolated = {{},std::vector<double>(numcoords,0)};
}


double expression::integrate(int physreg, int integrationorder) { return integrate(physreg, NULL, integrationorder); }
double expression::integrate(int physreg, expression meshdeform, int integrationorder) { return integrate(physreg, &meshdeform, integrationorder); }

double expression::integrate(int physreg, expression* meshdeform, int integrationorder)
{
    // Make sure this expression is scalar and the mesh
    // deformation expression has the right size.
    if (not(isscalar()))
    {
        std::cout << "Error in 'expression' object: cannot integrate a nonscalar expression" << std::endl;
        abort();
    }
    int problemdimension = universe::mymesh->getmeshdimension();
    if (meshdeform != NULL && (meshdeform->countcolumns() != 1 || meshdeform->countrows() < problemdimension))
    {
        std::cout << "Error in 'expression' object: mesh deformation expression has size " << meshdeform->countrows() << "x" << meshdeform->countcolumns() << " (expected " << problemdimension << "x1)" << std::endl;
        abort();
    }

    // Get only the disjoint regions with highest dimension elements:
    std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();

    // Multiharmonic expressions are not allowed.
    if (not(isharmonicone(selecteddisjregs)))
    {
        std::cout << "Error in 'expression' object: cannot integrate a multiharmonic expression (only constant harmonic 1)" << std::endl;
        abort();
    }
    if (meshdeform != NULL && not(meshdeform->isharmonicone(selecteddisjregs)))
    {
        std::cout << "Error in 'expression' object: the mesh deformation expression cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        abort();
    }
    
    universe::allowestimatorupdate(true);

    double integralvalue = 0;
    // Send the disjoint regions with same element type numbers together:
    disjointregionselector mydisjregselector(selecteddisjregs, {});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);

        // Get the Gauss points:
        int elementtypenumber = (universe::mymesh->getdisjointregions())->getelementtypenumber(mydisjregs[0]);
        gausspoints mygausspoints(elementtypenumber, integrationorder);
        std::vector<double> evaluationpoints = mygausspoints.getcoordinates();
        std::vector<double> weights = mygausspoints.getweights();

        // Loop on all total orientations (if required):
        bool isorientationdependent = isvalueorientationdependent(mydisjregs) || (meshdeform != NULL && meshdeform->isvalueorientationdependent(mydisjregs));
        elementselector myselector(mydisjregs, isorientationdependent);
        do
        {
            std::shared_ptr<jacobian> myjacobian(new jacobian(myselector, evaluationpoints, meshdeform));

            densematrix detjac = myjacobian->getdetjac();
            // The Jacobian determinant should be positive irrespective of the node numbering:
            detjac.abs();

            // Store it in the universe for reuse.
            universe::computedjacobian = myjacobian;
            universe::allowreuse();

            densematrix compxinterpolated = myoperations[0]->interpolate(myselector, evaluationpoints, meshdeform)[1][0];
            compxinterpolated.multiplyelementwise(detjac);

            universe::forbidreuse();

            densematrix weightsmat(mygausspoints.count(), 1, weights);
            compxinterpolated = compxinterpolated.multiply(weightsmat);

            integralvalue += compxinterpolated.sum();
        }
        while (myselector.next());
    }
    
    universe::allowestimatorupdate(false);
    
    return integralvalue;
}

void expression::write(int physreg, int numfftharms, std::string filename, int lagrangeorder)
{
    write(physreg, numfftharms, NULL, filename, lagrangeorder, -1);
}

void expression::write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder)
{
    write(physreg, numfftharms, &meshdeform, filename, lagrangeorder, -1);
}

void expression::write(int physreg, std::string filename, int lagrangeorder, int numtimesteps)
{
    write(physreg, -1, NULL, filename, lagrangeorder, numtimesteps);
}

void expression::write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps)
{
    write(physreg, -1, &meshdeform, filename, lagrangeorder, numtimesteps);
}


void expression::write(int physreg, int numfftharms, expression* meshdeform, std::string filename, int lagrangeorder, int numtimesteps)
{
    // Make sure this expression is a column vector and the
    // mesh deformation expression has the right size.
    if (mynumrows > 3 || mynumcols != 1)
    {
        std::cout << "Error in 'expression' object: can not write a " << mynumrows << "x" << mynumcols << " expression to file (only scalars or 2x1 or 3x1 column vectors)" << std::endl;
        abort();
    }
    int problemdimension = universe::mymesh->getmeshdimension();
    if (meshdeform != NULL && (meshdeform->countcolumns() != 1 || meshdeform->countrows() < problemdimension))
    {
        std::cout << "Error in 'expression' object: mesh deformation expression has size " << meshdeform->countrows() << "x" << meshdeform->countcolumns() << " (expected " << problemdimension << "x1)" << std::endl;
        abort();
    }

    // Minimum lagrange order is 1!
    if (lagrangeorder < 1)
        lagrangeorder = 1;

    field x("x"), y("y"), z("z");
    expression xyz(3,1, {x,y,z});

    // Loop on all disjoint regions:
    std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions();

    // Make sure the 'meshdeform' expression is constant in time:
    if (meshdeform != NULL && not(meshdeform->isharmonicone(selecteddisjregs)))
    {
        std::cout << "Error in 'expression' object: the mesh deformation expression cannot be multiharmonic (only constant harmonic 1)" << std::endl;
        abort();
    }

    // Get the geometry interpolation order (1 if the element is not curved):
    int geolagrangeorder = lagrangeorder;
    if (iointerface::isonlyisoparametric(filename) == false && meshdeform == NULL)
        geolagrangeorder = universe::mymesh->getelements()->getcurvatureorder();

    // These are the time tags that will be used:
    std::vector<double> timetags = {};
    if (numtimesteps > 0)
        timetags = std::vector<double>(numtimesteps, 0.0);

    // The data to write for harmonic h will be added at datatowrite[h][0].
    // For a time solution the data is at datatowrite[0][0].
    std::vector<std::vector<iodata>> datatowrite = {};
    if (numtimesteps > 0)
        datatowrite = {{iodata(lagrangeorder, geolagrangeorder, isscalar(), timetags)}};

    universe::allowestimatorupdate(true);

    // Send the disjoint regions with same element type numbers together:
    disjointregionselector mydisjregselector(selecteddisjregs, {});
    for (int g = 0; g < mydisjregselector.countgroups(); g++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(g);

        int elementtype = (universe::mymesh->getdisjointregions())->getelementtypenumber(mydisjregs[0]);
        element myelement(elementtype);

        // The expression will be interpolated at the following Lagrange nodes:
        lagrangeformfunction mylagrange(elementtype, lagrangeorder, {});
        std::vector<double> lagrangecoords = mylagrange.getnodecoordinates();
        // The x, y and z coordinates will be interpolated at the following Lagrange nodes:
        lagrangeformfunction mygeolagrange(elementtype, geolagrangeorder, {});
        std::vector<double> geolagrangecoords = mygeolagrange.getnodecoordinates();

        // Loop on all total orientations (if required):
        bool isorientationdependent = isvalueorientationdependent(mydisjregs) || (meshdeform != NULL && meshdeform->isvalueorientationdependent(mydisjregs));
        elementselector myselector(mydisjregs, isorientationdependent);
        do
        {
            // The harmonic content might have changed:
            std::vector<int> harms = {};

            // Compute the mesh coordinates. Initialise all coordinates to zero.
            std::vector<densematrix> coords(3,densematrix(myselector.countinselection(), geolagrangecoords.size()/3,0));
            for (int i = 0; i < problemdimension; i++)
            {
                coords[i] = (xyz.myoperations[i]->interpolate(myselector, geolagrangecoords, NULL))[1][0];
                if (meshdeform != NULL)
                    coords[i].add((meshdeform->myoperations[i]->interpolate(myselector, geolagrangecoords, NULL))[1][0]);
            }
            // Interpolate the current expression:
            std::vector<  std::vector<std::vector<densematrix>>  > expr(countrows());
            std::vector<densematrix> fftexpr(countrows());
            // Reuse what's possible to reuse during the interpolation:
            universe::allowreuse();
            for (int i = 0; i < countrows(); i++)
            {
                if (numtimesteps <= 0)
                {
                    if (numfftharms <= 0)
                        expr[i] = myoperations[i]->interpolate(myselector, lagrangecoords, meshdeform);
                    else
                        expr[i] = myfft::fft(myoperations[i]->multiharmonicinterpolate(numfftharms, myselector, lagrangecoords, meshdeform), myselector.countinselection(), lagrangecoords.size()/3);
                }
                else
                    fftexpr[i] = myfft::toelementrowformat(myoperations[i]->multiharmonicinterpolate(numtimesteps, myselector, lagrangecoords, meshdeform), myselector.countinselection());
            }
            universe::forbidreuse();

            // Make sure the harmonic content is the same for every component:
            if (numtimesteps <= 0)
                myfft::sameharmonics(expr);

            // Get a vector containing all harmonic numbers in 'expr' on the current disjoint regions:
            if (numtimesteps <= 0)
            {
                for (int h = 0; h < expr[0].size(); h++)
                {
                    if (expr[0][h].size() == 1)
                    {
                        harms.push_back(h);
                        if (datatowrite.size() < h+1)
                            datatowrite.resize(h+1);
                        if (datatowrite[h].size() == 0)
                            datatowrite[h] = {iodata(lagrangeorder, geolagrangeorder, isscalar(), timetags)};
                    }
                }
            }

            if (numtimesteps <= 0)
            {
                for (int h = 0; h < harms.size(); h++)
                {
                    datatowrite[harms[h]][0].addcoordinates(elementtype, coords[0], coords[1], coords[2]);
                    std::vector<densematrix> curdata(countrows());
                    for (int comp = 0; comp < countrows(); comp++)
                        curdata[comp] = expr[comp][harms[h]][0];
                    datatowrite[harms[h]][0].adddata(elementtype, curdata);
                }
            }
            else
            {
                datatowrite[0][0].addcoordinates(elementtype, coords[0], coords[1], coords[2]);
                datatowrite[0][0].adddata(elementtype, fftexpr);
            }
        }
        while (myselector.next());
    }
    
    universe::allowestimatorupdate(false);

    // Write the data:
    if (numtimesteps <= 0)
    {
        for (int h = 0; h < datatowrite.size(); h++)
        {
            if (datatowrite[h].size() == 1)
            {
                if (h <= 1)
                    iointerface::writetofile(filename, datatowrite[h][0]);
                else
                    iointerface::writetofile(filename, datatowrite[h][0], "_harm"+std::to_string(h));
            }
        }
    }
    else
        iointerface::writetofile(filename, datatowrite[0][0], "timesteps"+std::to_string(numtimesteps));
}

void expression::streamline(int physreg, std::string filename, const std::vector<double>& startcoords, double stepsize, bool downstreamonly)
{
    if (startcoords.size() == 0)
        return;

    // This can happen with int divisions:
    if (stepsize == 0)
    {
        std::cout << "Error in 'expression' object: step size for streamline cannot be zero" << std::endl;
        abort();
    }

    // Stream lines can only be obtained for expressions with at least as many components as the geometry dimension:
    int problemdimension = universe::mymesh->getmeshdimension();
    if (mynumrows < problemdimension || mynumrows > 3 || mynumcols != 1)
    {
        std::cout << "Error in 'expression' object: expected a column vector expression with " << problemdimension << " to 3 components to get the stream lines" << std::endl;
        abort();
    }
    // For simplicity the code below is written only for 3x1 expressions:
    if (mynumrows == 1)
    {
        expression(3,1,{this->getcopy(),0,0}).streamline(physreg, filename, startcoords, stepsize, downstreamonly);
        return;
    }
    if (mynumrows == 2)
    {
        expression(3,1,{this->at(0,0),this->at(1,0),0}).streamline(physreg, filename, startcoords, stepsize, downstreamonly);
        return;
    }

    if (startcoords.size()%3 != 0)
    {
        std::cout << "Error in 'expression' object: expected a vector with a length multiple of 3 for the stream line starting coordinates" << std::endl;
        abort();
    }

    int numnodes = startcoords.size()/3;

    // If upstream AND downstream calculation there are double the amount of starting coords:
    int factortwo = 2;
    if (downstreamonly)
        factortwo = 1;

    std::vector<std::vector<double>> xcoords(factortwo*numnodes), ycoords(factortwo*numnodes), zcoords(factortwo*numnodes), magnitude(factortwo*numnodes);

    // Current coordinates:
    std::vector<double> curcoords(factortwo* 3*numnodes);
    // Know which stream line is still active:
    std::vector<bool> isactive(factortwo*numnodes, true);

    // Stepsize vector for every starting point:
    std::vector<double> h(factortwo*numnodes);


    for (int i = 0; i < numnodes; i++)
    {
        curcoords[factortwo*3*i+0] = startcoords[3*i+0];
        curcoords[factortwo*3*i+1] = startcoords[3*i+1];
        curcoords[factortwo*3*i+2] = startcoords[3*i+2];

        h[factortwo*i+0] = stepsize;

        if (downstreamonly == false)
        {
            curcoords[factortwo*3*i+3] = startcoords[3*i+0];
            curcoords[factortwo*3*i+4] = startcoords[3*i+1];
            curcoords[factortwo*3*i+5] = startcoords[3*i+2];

            h[factortwo*i+1] = -stepsize;
        }
    }

    std::vector<double> k1, k2, k3, k4, k1scaled, k2scaled, k3scaled, k4scaled;
    std::vector<bool> isfound;

    bool isdone = false;
    while (isdone == false)
    {
        // Calculate the Runge-Kutta parameters k1, k2, k3 and k4:
        interpolate(physreg, NULL, curcoords, k1, isfound);
        k1scaled = myalgorithm::normblocks(k1,3);
        std::vector<double> ynplushk1over2 = curcoords;
        for (int i = 0; i < curcoords.size(); i++)
            ynplushk1over2[i] += h[(i-i%3)/3]*0.5*k1scaled[i];
        interpolate(physreg, NULL, ynplushk1over2, k2, isfound);
        k2scaled = myalgorithm::normblocks(k2,3);
        std::vector<double> ynplushk2over2 = curcoords;
        for (int i = 0; i < curcoords.size(); i++)
            ynplushk2over2[i] += h[(i-i%3)/3]*0.5*k2scaled[i];
        interpolate(physreg, NULL, ynplushk2over2, k3, isfound);
        k3scaled = myalgorithm::normblocks(k3,3);
        std::vector<double> ynplushk3 = curcoords;
        for (int i = 0; i < curcoords.size(); i++)
            ynplushk3[i] += h[(i-i%3)/3]*k3scaled[i];
        interpolate(physreg, NULL, ynplushk3, k4, isfound);
        k4scaled = myalgorithm::normblocks(k4,3);


        isdone = true;
        for (int i = 0; i < isfound.size(); i++)
        {
            if (isfound[i] == false)
                isactive[i] = false;
            if (isactive[i])
                isdone = false;
        }

        // Append data to write to disk:
        for (int i = 0; i < isactive.size(); i++)
        {
            if (isactive[i])
            {
                xcoords[i].push_back(curcoords[3*i+0]);
                ycoords[i].push_back(curcoords[3*i+1]);
                zcoords[i].push_back(curcoords[3*i+2]);

                double flowspeed = std::sqrt(k1[3*i+0]*k1[3*i+0] + k1[3*i+1]*k1[3*i+1] + k1[3*i+2]*k1[3*i+2]);

                magnitude[i].push_back(flowspeed);
            }
        }

        // Update the coordinates:
        for (int i = 0; i < curcoords.size(); i++)
            curcoords[i] += h[(i-i%3)/3]/6.0*( k1scaled[i]+2.0*k2scaled[i]+2.0*k3scaled[i]+k4scaled[i] );
    }


    // Write to file:
    iodata datatowrite(1, 1, true, {});
    for (int m = 0; m < xcoords.size(); m++)
    {
        int numlines = xcoords[m].size();
        if (numlines == 0)
            continue;

        densematrix xcoordsmat(numlines-1,2, 0), ycoordsmat(numlines-1,2, 0), zcoordsmat(numlines-1,2, 0), flowspeedmat(numlines-1,2, 0);
        double* xptr = xcoordsmat.getvalues(); double* yptr = ycoordsmat.getvalues(); double* zptr = zcoordsmat.getvalues();
        double* fsptr = flowspeedmat.getvalues();

        for (int i = 0; i < numlines-1; i++)
        {
            xptr[2*i+0] = xcoords[m][i%numlines]; xptr[2*i+1] = xcoords[m][(i+1)%numlines];
            yptr[2*i+0] = ycoords[m][i%numlines]; yptr[2*i+1] = ycoords[m][(i+1)%numlines];
            zptr[2*i+0] = zcoords[m][i%numlines]; zptr[2*i+1] = zcoords[m][(i+1)%numlines];

            fsptr[2*i+0] = magnitude[m][i%numlines]; fsptr[2*i+1] = magnitude[m][(i+1)%numlines];
        }

        datatowrite.addcoordinates(1, xcoordsmat, ycoordsmat, zcoordsmat);
        datatowrite.adddata(1, {flowspeedmat});
    }
    iointerface::writetofile(filename, datatowrite);
}

void expression::reuseit(bool istobereused)
{
    for (int i = 0; i < myoperations.size(); i++)
    {
        if (myoperations[i]->isdofincluded() == false && myoperations[i]->istfincluded() == false)
            myoperations[i]->reuseit(istobereused);
        else
        {
            std::cout << "Error in 'expression' object: cannot reuse an expression including a dof() or tf()" << std::endl;
            abort();
        }
    }
}

bool expression::isharmonicone(std::vector<int> disjregs)
{
    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (myoperations[i]->isharmonicone(disjregs) == false)
            return false;
    }
    return true;
}

bool expression::isvalueorientationdependent(std::vector<int> disjregs)
{
    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (myoperations[i]->isvalueorientationdependent(disjregs))
            return true;
    }
    return false;
}

bool expression::iszero(void)
{
    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (myoperations[i]->isconstant() == false || myoperations[i]->getvalue() != 0)
            return false;
    }
    return true;
}

vec expression::atbarycenter(int physreg, field onefield)
{
    // The field must be a "one" type:
    if (onefield.getpointer()->gettypename() != "one")
    {
        std::cout << "Error in 'expression' object: atbarycenter requires a 'one' type field" << std::endl;
        abort();
    }
    if (countcolumns() != 1 || onefield.countcomponents() != countrows())
    {
        std::cout << "Error in 'expression' object: in atbarycenter the size of the expression and of the argument field must match" << std::endl;
        abort();
    }

    // Compute the expression at the barycenter:
    formulation formul;
    
    integration myterm(physreg, - mathop::tf(onefield)*(this->getcopy()));
    myterm.isbarycentereval = true;

    formul += myterm;
    formul.generate();

    return formul.rhs();
}

void expression::print(void)
{
    std::cout << std::endl << "Expression size is " << mynumrows << "x" << mynumcols << std::endl;

    for (int row = 0; row < mynumrows; row++)
    {
        for (int col = 0; col < mynumcols; col++)
        {
            std::cout << " @ row " << row << ", col " << col << " :" << std::endl;
            myoperations[row*mynumcols+col]->print();
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

void expression::rotate(double ax, double ay, double az, std::string leftop, std::string rightop)
{
    if (isscalar())
    {
        std::cout << "Error in 'expression' object: cannot rotate a scalar expression" << std::endl;
        abort();
    }
	
    // Correctly transform a 3x3 and 3x1 by default:
    if (leftop == "default")
        leftop = "R";
    if (rightop == "default")
    {
        if (mynumcols > 1)
            rightop = "RT";
        else
            rightop = "";
    }
    

    if (leftop != "" && leftop != "R" && leftop != "RT" && leftop != "R-1" && leftop != "R-T" && leftop != "K" && leftop != "KT" && leftop != "K-1" && leftop != "K-T")
    {
        std::cout << "Error in 'expression' object: in 'rotate' left product can only be '', 'R', 'RT', 'R-1', 'R-T', 'K', 'KT', 'K-1', 'K-T'" << std::endl;
        abort();
    }
    if (rightop != "" && rightop != "R" && rightop != "RT" && rightop != "R-1" && rightop != "R-T" && rightop != "K" && rightop != "KT" && rightop != "K-1" && rightop != "K-T")
    {
        std::cout << "Error in 'expression' object: in 'rotate' right product can only be '', 'R', 'RT', 'R-1', 'R-T', 'K', 'KT', 'K-1', 'K-T'" << std::endl;
        abort();
    }
	
    // Define the rotation matrices needed:
    expression R,RT,K,KT,invK,invKT;
    
    if (leftop != "" && leftop[0] == 'R' || rightop != "" && rightop[0] == 'R')
    {
        R = mathop::rotation(ax, ay, az)[0];
        RT = mathop::transpose(R);   
    }
        
    if (leftop != "" && leftop[0] == 'K' || rightop != "" && rightop[0] == 'K')
    {
        std::vector<expression> Ks = mathop::rotation(ax, ay, az, "voigt");
        K = Ks[0]; invK = Ks[1];
        KT = mathop::transpose(K);
        invKT = mathop::transpose(invK);
    }
    
    
    ///// Rotate the matrix in this expression:
    
    expression rotated = this->getcopy();
    
    if (leftop == "R")
        rotated = R*rotated;
    if (leftop == "RT")
        rotated = RT*rotated;
    if (leftop == "R-1")
        rotated = RT*rotated;
    if (leftop == "R-T")
        rotated = R*rotated;
    if (leftop == "K")
        rotated = K*rotated;
    if (leftop == "KT")
        rotated = KT*rotated;
    if (leftop == "K-1")
        rotated = invK*rotated;
    if (leftop == "K-T")
        rotated = invKT*rotated;

    if (rightop == "R")
        rotated = rotated*R;
    if (rightop == "RT")
        rotated = rotated*RT;
    if (rightop == "R-1")
        rotated = rotated*RT;
    if (rightop == "R-T")
        rotated = rotated*R;
    if (rightop == "K")
        rotated = rotated*K;
    if (rightop == "KT")
        rotated = rotated*KT;
    if (rightop == "K-1")
        rotated = rotated*invK;
    if (rightop == "K-T")
        rotated = rotated*invKT;
        
    for (int i = 0; i < rotated.myoperations.size(); i++)
        rotated.myoperations[i] = rotated.myoperations[i]->simplify({});
    
    myoperations = rotated.myoperations;
}

expression expression::at(int row, int col)
{
    expression arrayentry;

    if (mynumrows < row+1 || mynumcols < col+1)
    {
        std::cout << "Error in 'expression' object: trying to get entry (" << row << ", " << col << ") in a " << mynumrows << " by " << mynumcols << " expression array" << std::endl;
        abort();
    }
    arrayentry.mynumrows = 1;
    arrayentry.mynumcols = 1;
    arrayentry.myoperations = {myoperations[row*mynumcols+col]};

    return arrayentry;
}

std::vector<double> expression::evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords)
{
    if (isscalar())
    {
        if (xcoords.size() == ycoords.size() && xcoords.size() == zcoords.size())
        {
            if (xcoords.size() == 0)
                return std::vector<double>(0);
            else
                return myoperations[0]->evaluate(xcoords, ycoords, zcoords);
        }
        else
        {
            std::cout << "Error in 'expression' object: expected vectors of same length as arguments in 'evaluate'" << std::endl;
            abort();
        }
    }
    else
    {
        std::cout << "Error in 'expression' object: 'evaluate' can only be called on a scalar expression" << std::endl;
        abort();
    }
}

expression expression::resize(int numrows, int numcols)
{
    std::vector<expression> args(numrows*numcols);
    for (int i = 0; i < numrows; i++)
    {
        for (int j = 0; j < numcols; j++)
        {
            if (i < mynumrows && j < mynumcols)
                args[i*numcols+j] = at(i,j);
            else
                args[i*numcols+j] = expression(0);
        }
    }

    expression output(numrows,numcols,args);
    return output;
}


expression expression::spacederivative(int whichderivative)
{
    if (mynumrows > 3 || mynumcols != 1)
    {
        std::cout << "Error in 'expression' object: can only take space derivatives of column vectors with up to 3 components" << std::endl;
        abort();
    }

    int problemdimension = universe::mymesh->getmeshdimension();

    expression derivated = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (derivated.myoperations[i]->isconstant())
            derivated.myoperations[i] = std::shared_ptr<operation>(new opconstant(0));
        else
        {
            if (whichderivative > problemdimension)
            {
                std::shared_ptr<opproduct> op(new opproduct( {derivated.myoperations[i], std::shared_ptr<operation>(new opconstant(0))} ));
                derivated.myoperations[i] = op;
            }
            else
            {
                std::shared_ptr<operation> op = derivated.myoperations[i]->copy();
                op->setspacederivative(whichderivative);
                derivated.myoperations[i] = op;
            }
        }
    }

    return derivated;
}

expression expression::kietaphiderivative(int whichderivative)
{
    expression derivated = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        std::shared_ptr<operation> op = derivated.myoperations[i]->copy();
        op->setkietaphiderivative(whichderivative);
        derivated.myoperations[i] = op;
    }

    return derivated;
}

expression expression::timederivative(int derivativeorder)
{
    expression derivated = this->getcopy();
    if (inrefcoord.size() > 0)
        derivated = inrefcoord[0].second;

    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (derivated.myoperations[i]->isconstant())
            derivated.myoperations[i] = std::shared_ptr<operation>(new opconstant(0));
        else
        {
            std::shared_ptr<operation> op = derivated.myoperations[i]->copy();
            op->increasetimederivativeorder(derivativeorder);
            derivated.myoperations[i] = op;
        }
    }

    // Transform from the reference element to the physical one a hcurl field:
    if (inrefcoord.size() > 0 && inrefcoord[0].first == "hcurl")
    {
        derivated.inrefcoord = {std::make_pair("hcurl",derivated)};
        derivated.myoperations = ( invjac()*derivated ).myoperations;
    }
    return derivated;
}

std::vector<std::pair<std::string,expression>> expression::getinrefcoord(void)
{
    return inrefcoord;
}

expression expression::transpose(void)
{
    expression transposed;

    // Inverted!
    transposed.mynumrows = mynumcols;
    transposed.mynumcols = mynumrows;

    transposed.myoperations.resize(mynumrows*mynumcols);

    for (int row = 0; row < transposed.mynumrows; row++)
    {
        for (int col = 0; col < transposed.mynumcols; col++)
            transposed.myoperations[row*transposed.mynumcols+col] = myoperations[col*mynumcols+row];
    }
    return transposed;
}

expression expression::removerowandcol(int rowtoremove, int coltoremove)
{
    expression submatrix;
    submatrix.mynumrows = mynumrows-1;
    submatrix.mynumcols = mynumcols-1;
    submatrix.myoperations.resize((mynumrows-1)*(mynumcols-1));

    int subrow = 0;
    for (int row = 0; row < mynumrows; row++)
    {
        if (row == rowtoremove) { continue; }
        int subcol = 0;
        for (int col = 0; col < mynumcols; col++)
        {
            if (col == coltoremove) { continue; }
            submatrix.myoperations[subrow*(mynumcols-1)+subcol] = myoperations[row*mynumcols+col];
            subcol++;
        }
        subrow++;
    }
    return submatrix;
}

expression expression::determinant(void)
{
    if (mynumrows == 1 && mynumcols == 1)
        return this->getcopy();

    if (mynumrows != mynumcols)
    {
        std::cout << "Error in 'expression' object: cannot get the determinant of a non-square matrix" << std::endl;
        abort();
    }

    std::shared_ptr<opsum> determ(new opsum());

    int row = 0;
    for (int col = 0; col < mynumcols; col++)
    {
        std::shared_ptr<operation> subdeterm = removerowandcol(row, col).determinant().myoperations[0];

        if ((row+col)%2 == 0)
            determ->addterm( std::shared_ptr<opproduct>(new opproduct( {myoperations[row*mynumcols+col], subdeterm} )) );
        else
            determ->addterm( std::shared_ptr<opproduct>(new opproduct( {std::shared_ptr<opconstant>(new opconstant(-1)), myoperations[row*mynumcols+col], subdeterm} )) );
    }

    expression output;
    output.mynumrows = 1;
    output.mynumcols = 1;
    output.myoperations = {determ};
    return output;
}

expression expression::cofactormatrix(void)
{
    expression cofactors = this->getcopy();

    for (int row = 0; row < mynumrows; row++)
    {
        for (int col = 0; col < mynumcols; col++)
        {
            std::shared_ptr<operation> op = removerowandcol(row,col).determinant().myoperations[0];

            if ((row+col)%2 == 0)
                cofactors.myoperations[row*mynumcols+col] = op;
            else
                cofactors.myoperations[row*mynumcols+col] = std::shared_ptr<opproduct>(new opproduct( {std::shared_ptr<opconstant>(new opconstant(-1)), op}));
        }
    }
    return cofactors;
}

expression expression::invert(void)
{
    if (mynumrows == 1 && mynumcols == 1)
        return 1/(this->getcopy());

    if (mynumrows != mynumcols)
    {
        std::cout << "Error in 'expression' object: cannot invert a non-square matrix" << std::endl;
        abort();
    }

    expression invdet = 1/determinant();
    // The determinant inverse will be reused a lot:
    invdet.reuseit();

    return invdet * cofactormatrix().transpose();
}

expression expression::pow(expression input)
{
    expression powered = this->getcopy();

    // The base and exponent expressions must be scalar:
    if (not(isscalar()) || not(input.isscalar()))
    {
        std::cout << "Error in 'expression' object: cannot have non scalar arguments in a power" << std::endl;
        abort();
    }

    std::shared_ptr<operation> myexponent = input.myoperations[0];
    std::shared_ptr<operation> mybase = myoperations[0];
    // In case the exponent is 0 or 1 we simplify.
    // Remember that doubles can exactly represent integers up to a large enough value.
    if (myexponent->isconstant() && (myexponent->getvalue() == 0 || myexponent->getvalue() == 1))
    {
        double val = myexponent->getvalue();
        // With a 0 power we get a constant. With a 1 power we get the base.
        if (val == 0)
            powered.myoperations = {std::shared_ptr<opconstant>(new opconstant(1))};
        else
            powered.myoperations = {mybase};
    }
    else
        powered.myoperations = {std::shared_ptr<oppower>(new oppower(mybase, myexponent))};

    return powered;
}

expression expression::dof(int physreg)
{
    expression doftag = this->getcopy();
    if (inrefcoord.size() > 0)
        doftag = inrefcoord[0].second;

    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (doftag.myoperations[i]->isconstant() && doftag.myoperations[i]->getvalue() == 0)
            continue;

        if (doftag.myoperations[i]->isfield())
        {
            if (doftag.myoperations[i]->getspacederivative() != 0 || doftag.myoperations[i]->getkietaphiderivative() != 0 || doftag.myoperations[i]->gettimederivative() != 0)
            {
                std::cout << "Error in 'expression' object: cannot apply space or time derivatives to the dof() field argument" << std::endl;
                abort();
            }
            std::shared_ptr<opdof> op(new opdof(doftag.myoperations[i]->getfieldpointer(), physreg));
            op->selectformfunctioncomponent(doftag.myoperations[i]->getformfunctioncomponent());
            op->setfieldcomponent(doftag.myoperations[i]->getfieldcomponent());
            doftag.myoperations[i] = op;
        }
        else
        {
            std::cout << "Error in 'expression' object: the argument of dof() must be a field or a field expression (constant 0 is allowed)" << std::endl;
            abort();
        }
    }

    // Transform from the reference element to the physical one a hcurl field:
    if (inrefcoord.size() > 0 && inrefcoord[0].first == "hcurl")
    {
        doftag.inrefcoord = {std::make_pair("hcurl",doftag)};
        doftag.myoperations = ( invjac()*doftag ).myoperations;
    }
    return doftag;
}

expression expression::tf(int physreg)
{
    expression tftag = this->getcopy();
    if (inrefcoord.size() > 0)
        tftag = inrefcoord[0].second;
        
    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (tftag.myoperations[i]->isconstant() && tftag.myoperations[i]->getvalue() == 0)
            continue;

        if (tftag.myoperations[i]->isfield())
        {
            if (tftag.myoperations[i]->getspacederivative() != 0 || tftag.myoperations[i]->getkietaphiderivative() != 0 || tftag.myoperations[i]->gettimederivative() != 0)
            {
                std::cout << "Error in 'expression' object: cannot apply space or time derivatives to the tf() field argument" << std::endl;
                abort();
            }
            std::shared_ptr<optf> op(new optf(tftag.myoperations[i]->getfieldpointer(), physreg));
            op->selectformfunctioncomponent(tftag.myoperations[i]->getformfunctioncomponent());
            op->setfieldcomponent(tftag.myoperations[i]->getfieldcomponent());
            tftag.myoperations[i] = op;
        }
        else
        {
            std::cout << "Error in 'expression' object: the argument of tf() must be a field or a field expression (constant 0 is allowed)" << std::endl;
            abort();
        }
    }

    // Transform from the reference element to the physical one a hcurl field:
    if (inrefcoord.size() > 0 && inrefcoord[0].first == "hcurl")
    {
        tftag.inrefcoord = {std::make_pair("hcurl",tftag)};
        tftag.myoperations = ( invjac()*tftag ).myoperations;
    }
    return tftag;
}

expression expression::sin(void)
{
    expression sinexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        sinexpr.myoperations[i] = std::shared_ptr<opsin>(new opsin(myoperations[i]));

    return sinexpr;
}

expression expression::cos(void)
{
    expression cosexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        cosexpr.myoperations[i] = std::shared_ptr<opcos>(new opcos(myoperations[i]));

    return cosexpr;
}

expression expression::tan(void)
{
    expression tanexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        tanexpr.myoperations[i] = std::shared_ptr<optan>(new optan(myoperations[i]));

    return tanexpr;
}

expression expression::asin(void)
{
    expression asinexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        asinexpr.myoperations[i] = std::shared_ptr<opasin>(new opasin(myoperations[i]));

    return asinexpr;
}

expression expression::acos(void)
{
    expression acosexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        acosexpr.myoperations[i] = std::shared_ptr<opacos>(new opacos(myoperations[i]));

    return acosexpr;
}

expression expression::atan(void)
{
    expression atanexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        atanexpr.myoperations[i] = std::shared_ptr<opatan>(new opatan(myoperations[i]));

    return atanexpr;
}

expression expression::abs(void)
{
    expression absexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        absexpr.myoperations[i] = std::shared_ptr<opabs>(new opabs(myoperations[i]));

    return absexpr;
}

expression expression::log10(void)
{
    expression log10expr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        log10expr.myoperations[i] = std::shared_ptr<oplog10>(new oplog10(myoperations[i]));

    return log10expr;
}

expression expression::mod(double modval)
{
    expression modexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
        modexpr.myoperations[i] = std::shared_ptr<opmod>(new opmod(myoperations[i], modval));

    return modexpr;
}

expression expression::on(int physreg, expression* coordshift, bool errorifnotfound)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    if (coordshift != NULL && (coordshift->countcolumns() != 1 || coordshift->countrows() < problemdimension))
    {
        std::cout << "Error in 'expression' object: coordinate shift argument in 'on' has size " << coordshift->countrows() << "x" << coordshift->countcolumns() << " (expected " << problemdimension << "x1)" << std::endl;
        abort();
    }

    expression onexpr = this->getcopy();

    for (int i = 0; i < mynumrows*mynumcols; i++)
    {
        if (myoperations[i]->istfincluded())
        {
            std::cout << "Error in 'expression' object: argument of 'on' cannot include a test function tf()" << std::endl;
            abort();
        }
        if (myoperations[i]->isdofincluded() == false)
            onexpr.myoperations[i] = std::shared_ptr<opon>(new opon(physreg, coordshift, onexpr.myoperations[i], errorifnotfound));
        else
        {
            // Isolate the dofs (multiply by a dummy scalar test function for the call to 'extractdoftfpolynomial').
            field dummy("h1");
            expression curexpr = expression(onexpr.myoperations[i]) * mathop::tf(dummy);
            curexpr.expand();
        
            int elementdimension = universe::mymesh->getphysicalregions()->get(physreg)->getelementdimension();
            std::vector< std::vector<std::vector<std::shared_ptr<operation>>> > coeffdoftf = curexpr.extractdoftfpolynomial(elementdimension);
            // Do not retrieve the info for the tf:
            std::vector<std::vector<std::shared_ptr<operation>>> coeffs = coeffdoftf[0]; 
            std::vector<std::vector<std::shared_ptr<operation>>> dofs = coeffdoftf[1];
            
            std::vector<std::shared_ptr<operation>> allterms = {};
            for (int m = 0; m < dofs.size(); m++)
            {
                for (int n = 0; n < dofs[m].size(); n++)
                {
                    // The coefficient is a new opon object:
                    std::shared_ptr<operation> curcoef(new opon(physreg, coordshift, coeffs[m][n]->copy(), errorifnotfound));
                    // The dof gets the on tag:
                    std::shared_ptr<operation> curdof, curterm;
                    if (dofs[m][n]->getfieldpointer() != NULL)
                    {
                        curdof = dofs[m][n]->copy();
                        oncontext ctxt(physreg, coordshift, errorifnotfound);
                        curdof->setoncontext(ctxt);
                        curterm = std::shared_ptr<opproduct>(new opproduct({curcoef,curdof}));
                    }
                    else
                        curterm = curcoef;
			
                    curterm = curterm->simplify({});
			
                    allterms.push_back(curterm);
                }
            }
            if (allterms.size() == 0)
                allterms = {std::shared_ptr<opconstant>(new opconstant(0.0))};
            onexpr.myoperations[i] = std::shared_ptr<opsum>(new opsum(allterms));
        }
    }

    return onexpr;
}

expression expression::time(void)
{
    expression exp;
    exp.mynumrows = 1;
    exp.mynumcols = 1;
    exp.myoperations = {std::shared_ptr<optime>(new optime)};

    return exp;
}

expression expression::invjac(int row, int col)
{
    expression exp;
    exp.mynumrows = 1;
    exp.mynumcols = 1;
    exp.myoperations = {std::shared_ptr<opinvjac>(new opinvjac(row,col))};

    return exp;
}

expression expression::jac(int row, int col)
{
    expression exp;
    exp.mynumrows = 1;
    exp.mynumcols = 1;
    exp.myoperations = {std::shared_ptr<opjac>(new opjac(row,col))};

    return exp;
}

expression expression::detjac(void)
{
    expression exp;
    exp.mynumrows = 1;
    exp.mynumcols = 1;
    exp.myoperations = {std::shared_ptr<opdetjac>(new opdetjac)};

    return exp;
}

expression expression::invjac(void)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    switch (problemdimension)
    {
        case 1:
            return expression(3,3,{invjac(0,0),0,0,   0,1,0,   0,0,1});
        case 2:
        {
            if (universe::isaxisymmetric)
                return expression(3,3,{invjac(0,0),invjac(0,1),0,   invjac(1,0),invjac(1,1),0,   0,0,invjac(2,2)});
            else
                return expression(3,3,{invjac(0,0),invjac(0,1),0,   invjac(1,0),invjac(1,1),0,   0,0,1});
        }
        case 3:
            return expression(3,3,{invjac(0,0),invjac(0,1),invjac(0,2),   invjac(1,0),invjac(1,1),invjac(1,2),   invjac(2,0),invjac(2,1),invjac(2,2)});
    }
}

expression expression::jac(void)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    switch (problemdimension)
    {
        case 1:
            return expression(3,3,{jac(0,0),0,0,   0,1,0,   0,0,1});
        case 2:
        {
            if (universe::isaxisymmetric)
                return expression(3,3,{jac(0,0),jac(0,1),0,   jac(1,0),jac(1,1),0,   0,0,jac(2,2)});
            else
                return expression(3,3,{jac(0,0),jac(0,1),0,   jac(1,0),jac(1,1),0,   0,0,1});
        }
        case 3:
            return expression(3,3,{jac(0,0),jac(0,1),jac(0,2),   jac(1,0),jac(1,1),jac(1,2),   jac(2,0),jac(2,1),jac(2,2)});
    }
}

expression expression::getcopy(void)
{
    expression output = *this;
    output.inrefcoord = {};

    return output;
}

std::shared_ptr<operation> expression::getoperationinarray(int row, int col)
{
    if (mynumrows < row+1 || mynumcols < col+1)
    {
        std::cout << "Error in 'expression' object: trying to get entry (" << row << ", " << col << ") in a " << mynumrows << " by " << mynumcols << " expression array" << std::endl;
        abort();
    }
    return myoperations[row*mynumcols+col];
}

void expression::expand(void)
{
    if (isscalar())
        myoperations[0] = myoperations[0]->expand();
    else
    {
        std::cout << "Error in 'expression' object: expand is only defined for scalar expressions" << std::endl;
        std::cout << "Did you try to define a non scalarisable formulation?" << std::endl;
        abort();
    }
}

std::vector< std::vector<std::vector<std::shared_ptr<operation>>> > expression::extractdoftfpolynomial(int elementdimension)
{
    // Simplify the operation:
    myoperations[0] = myoperations[0]->simplify({});
    
    if (myoperations[0]->isconstant() && myoperations[0]->getvalue() == 0)
        return std::vector< std::vector<std::vector<std::shared_ptr<operation>>> >(3, std::vector<std::vector<std::shared_ptr<operation>>>(0));

    // To make sure the format of the expression is ok.
    // Otherwise an error is thrown at the end.
    bool isformatok = true;

    // Every slice in coeffs, dofs and tfs corresponds to a unique
    // dof field-tf field combination (with a unique combination
    // of time derivative operator and form function component).
    std::vector<std::vector<std::shared_ptr<operation>>> coeffs = {};
    std::vector<std::vector<std::shared_ptr<operation>>> dofs = {};
    std::vector<std::vector<std::shared_ptr<operation>>> tfs = {};

    // Extract the sum terms from the operation:
    std::vector<std::shared_ptr<operation>> sumterms;
    if (myoperations[0]->issum())
        sumterms = myoperations[0]->getarguments();
    else
        sumterms = {myoperations[0]};

    // Loop on all elementary sum terms in the formulation:
    for (int i = 0; i < sumterms.size(); i++)
    {
        // Deal with the very specific case of a tf() term without
        // coefficient by multiplying it by a constant 1 to get a product:
        if (sumterms[i]->istf())
        {
            std::shared_ptr<opproduct> op(new opproduct);
            op->multiplybyterm(sumterms[i]);
            op->multiplybyterm(std::shared_ptr<operation>(new opconstant(1)));
            sumterms[i] = op;
        }
        // In a valid formulation term sumterms[i] must now always be a product:
        if (sumterms[i]->isproduct())
        {
            std::shared_ptr<operation> currentdof(new opdof(NULL));
            std::shared_ptr<operation> currenttf(new optf(NULL));

            // The coef is what remains after dof() and tf() removal.
            // We thus remove the dof and tf terms.
            std::shared_ptr<operation> currentcoef = sumterms[i]->copy();

            std::vector<std::shared_ptr<operation>> productterms = sumterms[i]->getarguments();

            for (int j = productterms.size()-1; j >= 0; j--)
            {
                // Remove the dof and tf terms to get the coefficient:
                if (productterms[j]->isdof())
                {
                    currentdof = productterms[j];
                    currentcoef->removeterm(j);
                }
                if (productterms[j]->istf())
                {
                    currenttf = productterms[j];
                    currentcoef->removeterm(j);
                }
            }

            // A coef without factors actually has value 1:
            if (currentcoef->count() == 0)
                currentcoef = std::shared_ptr<opconstant>(new opconstant(1));

            // Do some error checking.
            // Make sure there is no dof or tf in the coef.
            if (currentcoef->isdofincluded() || currentcoef->istfincluded())
                isformatok = false;
            if (currenttf->getfieldpointer() == NULL)
            {
                std::cout << "Error in 'expression' object: malformed expression provided to the formulation (the test function is missing)" << std::endl;
                std::cout << "Expression was:" << std::endl;
                myoperations[0]->print();
                std::cout << std::endl;
                abort();
            }

            // Know which slice (i.e. dof field-tf field combination) we are at.
            // In a given slice all dofs and tfs must be defined on the same physical region!
            int currentslice = -1;
            for (int slice = 0; slice < tfs.size(); slice++)
            {
                // The pointed field, the physical region and the time derivative must be identical:
                if (dofs[slice][0]->getfieldpointer() == currentdof->getfieldpointer() && tfs[slice][0]->getfieldpointer() == currenttf->getfieldpointer() && dofs[slice][0]->getphysicalregion() == currentdof->getphysicalregion() && tfs[slice][0]->getphysicalregion() == currenttf->getphysicalregion() && dofs[slice][0]->gettimederivative() == currentdof->gettimederivative() && tfs[slice][0]->gettimederivative() == currenttf->gettimederivative())
                {
                    // Make sure the on context is the same (if any):
                    if (dofs[slice][0]->getfieldpointer() == NULL || ( dofs[slice][0]->getoncontext()->isequal(currentdof->getoncontext()) ))
                    {
                        currentslice = slice;
                        break;
                    }
                }
            }
            // If the slice does not yet exist create it.
            if (currentslice == -1)
            {
                currentslice = tfs.size();

                coeffs.push_back({});
                dofs.push_back({});
                tfs.push_back({});
            }

            // Explode any x, y, z space derivatived dof and tf into a sum of
            // products of invjac terms and ki, eta and phi derivatives.
            int numdofterms = 1, numtfterms = 1;
            if (currentdof->getspacederivative() > 0)
                numdofterms = elementdimension;
            if (currenttf->getspacederivative() > 0)
                numtfterms = elementdimension;

            std::vector<std::shared_ptr<operation>> currentdofsplit(numdofterms*numtfterms);
            std::vector<std::shared_ptr<operation>> currenttfsplit(numdofterms*numtfterms);
            std::vector<std::shared_ptr<operation>> currentcoefsplit(numdofterms*numtfterms);
            // Loop on every term in the product of the two sums:
            for (int i = 0; i < numdofterms; i++)
            {
                for (int j = 0; j < numtfterms; j++)
                {
                    std::shared_ptr<operation> copieddof = currentdof->copy();
                    std::shared_ptr<operation> copiedtf = currenttf->copy();
                    std::shared_ptr<operation> newcoef = currentcoef;
                    // Precompute the current coef when generating since it might appear several times:
                    currentcoef->reuseit(true);
                    // Multiply the coef by the appropriate invjac term:
                    if (currentdof->getspacederivative() > 0)
                    {
                        std::shared_ptr<opinvjac> invjacob(new opinvjac(currentdof->getspacederivative()-1, i));
                        newcoef = std::shared_ptr<operation>(new opproduct({newcoef, invjacob}));
                        copieddof->setkietaphiderivative(i+1);
                    }
                    if (currenttf->getspacederivative() > 0)
                    {
                        std::shared_ptr<opinvjac> invjacob(new opinvjac(currenttf->getspacederivative()-1, j));
                        newcoef = std::shared_ptr<operation>(new opproduct({newcoef, invjacob}));
                        copiedtf->setkietaphiderivative(j+1);
                    }
                    currentdofsplit[i*numtfterms+j] = copieddof;
                    currenttfsplit[i*numtfterms+j] = copiedtf;
                    currentcoefsplit[i*numtfterms+j] = newcoef;
                }
            }

            // Find the appropriate indexes in the current slice to put the split terms:
            for (int term = 0; term < currenttfsplit.size(); term++)
            {
                // In the current slice find the entry at which to
                // add the current coef, tf and dof (if existing).
                int currententry = -1;
                for (int entry = 0; entry < tfs[currentslice].size(); entry++)
                {
                    if (dofs[currentslice][entry]->getformfunctioncomponent() == currentdofsplit[term]->getformfunctioncomponent() && dofs[currentslice][entry]->getkietaphiderivative() == currentdofsplit[term]->getkietaphiderivative() && tfs[currentslice][entry]->getformfunctioncomponent() == currenttfsplit[term]->getformfunctioncomponent() && tfs[currentslice][entry]->getkietaphiderivative() == currenttfsplit[term]->getkietaphiderivative())
                    {
                        currententry = entry;
                        // Add the current coefficient to the entry:
                        std::shared_ptr<opsum> op(new opsum);
                        op->addterm(coeffs[currentslice][currententry]);
                        op->addterm(currentcoefsplit[term]);
                        coeffs[currentslice][currententry] = op;
                        break;
                    }
                }
                // If the entry does not yet exist create it.
                if (currententry == -1)
                {
                    currententry = tfs[currentslice].size();

                    coeffs[currentslice].push_back(currentcoefsplit[term]);
                    dofs[currentslice].push_back(currentdofsplit[term]);
                    tfs[currentslice].push_back(currenttfsplit[term]);
                }
            }
        }
        else
            isformatok = false;
    }

    if (not(isformatok))
    {
        std::cout << "Error in 'expression' object: don't know what to do with the expression provided to the formulation" << std::endl;
        std::cout << "The expression should be rewritable into a sum of products of the form coef*dof*tf (derivatives allowed)" << std::endl;
        std::cout << "Expression was:" << std::endl;
        myoperations[0]->print();
        std::cout << std::endl;
        abort();
    }

    return {coeffs, dofs, tfs};
}


expression expression::operator+(void) { return this->getcopy(); }
expression expression::operator-(void) { return (this->getcopy() * -1.0); }

expression expression::operator+(expression input)
{
    expression output = this->getcopy();

    if (mynumrows != input.mynumrows || mynumcols != input.mynumcols)
    {
        std::cout << "Error in 'expression' object: trying to add a " << mynumrows << "x" << mynumcols << " expression to a " << input.mynumrows << "x" << input.mynumcols << std::endl;
        abort();
    }
    for (int i = 0; i < mynumrows*mynumcols; i++)
        output.myoperations[i] = std::shared_ptr<opsum>(new opsum( {myoperations[i], input.myoperations[i]} ));
    return output;
}

expression expression::operator-(expression input) { return (this->getcopy() + input*-1); }

expression expression::operator*(expression input)
{
    expression output;

    // Products by scalars are always allowed:
    if (isscalar())
    {
        output = input;

        for (int i = 0; i < input.mynumrows*input.mynumcols; i++)
            output.myoperations[i] = std::shared_ptr<opproduct>(new opproduct( {myoperations[0], input.myoperations[i]} ));
        return output;
    }
    if (input.isscalar())
    {
        output = this->getcopy();

        for (int i = 0; i < mynumrows*mynumcols; i++)
            output.myoperations[i] = std::shared_ptr<opproduct>(new opproduct( {myoperations[i], input.myoperations[0]} ));
        return output;
    }

    // Scalar products of two row or two column vectors are allowed:
    if ( (mynumrows == 1 || mynumcols == 1) && mynumrows == input.mynumrows && mynumcols == input.mynumcols)
    {
        std::shared_ptr<opsum> op(new opsum);
        for (int i = 0; i < std::max(mynumrows,mynumcols); i++)
            op->addterm( std::shared_ptr<opproduct>(new opproduct( {myoperations[i], input.myoperations[i]} )) );
        output.mynumrows = 1;
        output.mynumcols = 1;
        output.myoperations = {op};
        return output;
    }

    // For all other products sizes must match.
    if (mynumcols != input.mynumrows)
    {
        std::cout << "Error in 'expression' object: trying to multiply a " << mynumrows << "x" << mynumcols << " expression by a " << input.mynumrows << "x" << input.mynumcols << std::endl;
        abort();
    }

    output.mynumrows = mynumrows;
    output.mynumcols = input.mynumcols;
    output.myoperations.resize(mynumrows*input.mynumcols);

    for (int row = 0; row < mynumrows; row++)
    {
        for (int col = 0; col < input.mynumcols; col++)
        {
            // Build every element in the output matrix:
            std::shared_ptr<opsum> op(new opsum);

            for (int k = 0; k < mynumcols; k++)
                op->addterm( std::shared_ptr<opproduct>(new opproduct( {myoperations[row*mynumcols+k], input.myoperations[k*input.mynumcols+col]} )) );

            output.myoperations[row*output.mynumcols+col] = op;
        }
    }
    return output;
}

expression expression::operator/(expression input)
{
    expression output;

    // Only divisions by scalars are allowed:
    if (input.isscalar())
    {
        output = this->getcopy();

        for (int i = 0; i < mynumrows*mynumcols; i++)
        {
            std::shared_ptr<opinversion> opdiv(new opinversion(input.myoperations[0]));
            std::shared_ptr<opproduct> opprod(new opproduct({myoperations[i],opdiv}));
            output.myoperations[i] = opprod;
        }
        return output;
    }
    else
    {
        std::cout << "Error in 'expression' object: only divisions by scalars are allowed" << std::endl;
        abort();
    }
}

expression expression::operator+(field inputfield) { return this->getcopy() + (expression)inputfield; }
expression expression::operator-(field inputfield) { return this->getcopy() - (expression)inputfield; }
expression expression::operator*(field inputfield) { return this->getcopy() * (expression)inputfield; }
expression expression::operator/(field inputfield) { return this->getcopy() / (expression)inputfield; }

expression expression::operator+(double val) { return this->getcopy() + (expression)val; }
expression expression::operator-(double val) { return this->getcopy() - (expression)val; }
expression expression::operator*(double val) { return this->getcopy() * (expression)val; }
expression expression::operator/(double val) { return this->getcopy() / (expression)val; }

expression expression::operator+(parameter& param) { return this->getcopy() + (expression)param; }
expression expression::operator-(parameter& param) { return this->getcopy() - (expression)param; }
expression expression::operator*(parameter& param) { return this->getcopy() * (expression)param; }
expression expression::operator/(parameter& param) { return this->getcopy() / (expression)param; }


expression operator+(double val, expression expr) { return expr+val; }
expression operator-(double val, expression expr) { return -expr+val; }
expression operator*(double val, expression expr) { return expr*val; }
expression operator/(double val, expression expr) { return ( (expression)val ) / expr; }

expression operator+(field inputfield, expression expr) { return expr+inputfield; }
expression operator-(field inputfield, expression expr) { return -expr+inputfield; }
expression operator*(field inputfield, expression expr) { return expr*inputfield; }
expression operator/(field inputfield, expression expr) { return ( (expression)inputfield ) / expr; }

expression operator+(parameter& param, expression expr) { return expr+param; }
expression operator-(parameter& param, expression expr) { return -expr+param; }
expression operator*(parameter& param, expression expr) { return expr*param; }
expression operator/(parameter& param, expression expr) { return ( (expression)param ) / expr; }
