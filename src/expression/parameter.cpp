#include "parameter.h"


void parameter::synchronize(void)
{
    if (issynchronizing || universe::mymesh->getmeshnumber() == mymeshnumber)
        return;
    issynchronizing = true;    


    // Flush the structure:
    myoperations = std::vector<std::vector<std::shared_ptr<operation>>>(universe::mymesh->getdisjointregions()->count(), std::vector<std::shared_ptr<operation>>(1, NULL));
    maxopnum = -1;
    opnums = std::vector<int>(universe::mymesh->getdisjointregions()->count(),-1);

    // Rebuild the structure:
    for (int i = 0; i < mystructuretracker.size(); i++)
        set(mystructuretracker[i].first, mystructuretracker[i].second);
    
    
    mymeshnumber = universe::mymesh->getmeshnumber();
    issynchronizing = false;
}

void parameter::errorifundefined(std::vector<int> disjregs)
{
    synchronize();

    for (int i = 0; i < disjregs.size(); i++)
    {
        if (myoperations[disjregs[i]].size() == 1 && myoperations[disjregs[i]][0] == NULL)
        {
            std::cout << "Error in 'parameter' object: the parameter has not been defined on the requested region" << std::endl;
            abort();
        }
    }
}

std::vector<int> parameter::getopnums(std::vector<int> disjregs)
{
    synchronize();
    
    std::vector<int> output(disjregs.size());
    for (int i = 0; i < disjregs.size(); i++)
        output[i] = opnums[disjregs[i]];
    return output;
}

parameter::parameter(void) : myoperations((universe::mymesh->getdisjointregions())->count(), std::vector<std::shared_ptr<operation>>(1, NULL)), opnums((universe::mymesh->getdisjointregions())->count(),-1)
{ 
    mynumrows = 1; 
    mynumcols = 1; 
}

parameter::parameter(int numrows, int numcols) : myoperations((universe::mymesh->getdisjointregions())->count(), std::vector<std::shared_ptr<operation>>(numrows*numcols, NULL)), opnums((universe::mymesh->getdisjointregions())->count(),-1)
{ 
    mynumrows = numrows; 
    mynumcols = numcols; 
}

void parameter::set(int physreg, expression input)
{
    synchronize();
    
    // Keep track of the calls to 'set':
    if (issynchronizing == false)
        mystructuretracker.push_back(std::make_pair(physreg, input));
        
    if (mynumrows != input.countrows() || mynumcols != input.countcolumns())
    {
        std::cout << "Error in 'parameter' object: trying to set the " << mynumrows << "x" << mynumcols << " sized parameter to a size " << input.countrows() << "x" << input.countcolumns() << std::endl;
        abort();
    }
    
    maxopnum++;
    
    // Consider ALL disjoint regions in the physical region with (-1):
    std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);

    for (int i = 0; i < selecteddisjregs.size(); i++)
    {
        opnums[selecteddisjregs[i]] = maxopnum;
        for (int row = 0; row < mynumrows; row++)
        {
            for (int col = 0; col < mynumcols; col++)
            {
                std::shared_ptr<operation> op = input.getoperationinarray(row, col);
                // Make sure there is no dof or tf in the operation:
                if (op->isdofincluded() || op->istfincluded())
                {
                    std::cout << "Error in 'parameter' object: cannot set an expression containing a dof or a tf" << std::endl;
                    abort();
                }
                myoperations[selecteddisjregs[i]][row*mynumcols+col] = op;
            }
        }
    }
}

std::shared_ptr<operation> parameter::get(int disjreg, int row, int col)
{
    synchronize();
    
    errorifundefined({disjreg});
    return myoperations[disjreg][row*mynumcols+col];
}

int parameter::countrows(void)
{
    synchronize();
    
    return mynumrows;
}

int parameter::countcolumns(void)
{
    synchronize();
    
    return mynumcols;
}

parameterselectedregion parameter::operator|(int physreg)
{
    synchronize();
    
    return parameterselectedregion(this, physreg);
}

std::vector<std::vector<densematrix>> parameter::interpolate(int row, int col, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    synchronize();
    
    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();
    
    // Make sure the parameter has been defined:
    errorifundefined(alldisjregs);
    
    // We have stored in the universe data defined on a set of disjoint 
    // regions. Here we will interpolate on subsets of disjoint regions 
    // so the data stored cannot be used during these interpolations.
    // Since we do not want to discard the computed data we simply
    // set 'isreuseallowed' temporarily to false.
    bool wasreuseallowed = universe::isreuseallowed;
    universe::isreuseallowed = false;
    
    std::vector<std::vector<densematrix>> out = {};
    
    // Group disj. regs. with same operation number (and same element type number).
    disjointregionselector mydisjregselector(alldisjregs, {getopnums(alldisjregs)});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
    
        elemselect.selectdisjointregions(mydisjregs);
        if (elemselect.countinselection() == 0)
            continue;
        elementselector myselection = elemselect.extractselection();
        
        // IMPORTANT: Harmonic numbers can be different from one disj. reg. to the other:
        std::vector<std::vector<densematrix>> currentinterp = myoperations[mydisjregs[0]][row*mynumcols+col]->interpolate(myselection, evaluationcoordinates, meshdeform);
        
        // Preallocate the harmonics not yet in 'out':
        if (out.size() < currentinterp.size())
            out.resize(currentinterp.size());
        for (int h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1 && out[h].size() == 0)
                out[h] = {densematrix(elemselect.countincurrentorientation(), evaluationcoordinates.size()/3, 0)};
        }
        
        // Insert 'currentinterp' in 'out':
        for (int h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1)
                out[h][0].insertatrows(elemselect.getelementindexes(), currentinterp[h][0]);
        }
    }
    // Unselect the disjoint regions:
    elemselect.selectdisjointregions({});
    
    // Reset the reuse right:
    universe::isreuseallowed = wasreuseallowed;
    
    return out;
}

densematrix parameter::multiharmonicinterpolate(int row, int col, int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    synchronize();
    
    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();
    
    // Make sure the parameter has been defined:
    errorifundefined(alldisjregs);
    
    // We have stored in the universe data defined on a set of disjoint 
    // regions. Here we will interpolate on subsets of disjoint regions 
    // so the data stored cannot be used during these interpolations.
    // Since we do not want to discard the computed data we simply
    // set 'isreuseallowed' temporarily to false.
    bool wasreuseallowed = universe::isreuseallowed;
    universe::isreuseallowed = false;
    
    // Preallocate the output matrix:
    densematrix out(numtimeevals, elemselect.countincurrentorientation() * evaluationcoordinates.size()/3);
    
    // Group disj. regs. with same operation number (and same element type number).
    disjointregionselector mydisjregselector(alldisjregs, {getopnums(alldisjregs)});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
        
        elemselect.selectdisjointregions(mydisjregs);
        if (elemselect.countinselection() == 0)
            continue;
        elementselector myselection = elemselect.extractselection();
        
        std::vector<int> selectedelementindexes = elemselect.getelementindexes();
        std::vector<int> selectedcolumns(selectedelementindexes.size()*evaluationcoordinates.size()/3);
        for (int j = 0; j < selectedelementindexes.size(); j++)
        {
            for (int k = 0; k < evaluationcoordinates.size()/3; k++)
                selectedcolumns[j*evaluationcoordinates.size()/3+k] = selectedelementindexes[j]*evaluationcoordinates.size()/3+k;
        }
        out.insertatcolumns(selectedcolumns, myoperations[mydisjregs[0]][row*mynumcols+col]->multiharmonicinterpolate(numtimeevals, myselection, evaluationcoordinates, meshdeform));
    }
    // Unselect the disjoint regions:
    elemselect.selectdisjointregions({});
    
    // Reset the reuse right:
    universe::isreuseallowed = wasreuseallowed;
    
    return out;
}    
    
void parameter::simplify(int row, int col, int disjreg)
{
    synchronize();
    
    // Make sure the parameter has been defined:
    errorifundefined({disjreg});
    
    myoperations[disjreg][row*mynumcols+col] = myoperations[disjreg][row*mynumcols+col]->simplify({disjreg});
}

void parameter::print(void)
{
    synchronize();
    
    std::cout << std::endl;
    std::cout << "Printing parameter of size " << mynumrows << "x" << mynumcols;
    std::cout << std::endl;
    for (int i = 0; i < myoperations.size(); i++)
    {   
        if (myoperations[i][0] != NULL)
        {
            std::cout << "------------ Expression on disjoint region " << i << " ------------" << std::endl << std::endl;
            for (int row = 0; row < mynumrows; row++)
            {          
                for (int col = 0; col < mynumcols; col++)
                {      
                    std::cout << " @ row " << row << ", col " << col << " :" << std::endl;
                    myoperations[i][row*mynumcols+col]->print();
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }
}





vec parameter::atbarycenter(int physreg, field onefield)
{ return ((expression)*this).atbarycenter(physreg, onefield); }

std::vector<double> parameter::max(int physreg, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).max(physreg, refinement, xyzrange); }
std::vector<double> parameter::max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).max(physreg, meshdeform, refinement, xyzrange); }
std::vector<double> parameter::min(int physreg, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).min(physreg, refinement, xyzrange); }
std::vector<double> parameter::min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange)
{ return ((expression)*this).min(physreg, meshdeform, refinement, xyzrange); }

void parameter::interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{ ((expression)*this).interpolate(physreg, xyzcoord, interpolated, isfound); }
void parameter::interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound)
{ ((expression)*this).interpolate(physreg, meshdeform, xyzcoord, interpolated, isfound); }

std::vector<double> parameter::interpolate(int physreg, const std::vector<double> xyzcoord)
{ return ((expression)*this).interpolate(physreg, xyzcoord); }
std::vector<double> parameter::interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord)
{ return ((expression)*this).interpolate(physreg, meshdeform, xyzcoord); }
        
double parameter::integrate(int physreg, expression meshdeform, int integrationorder) { return ((expression)*this).integrate(physreg, meshdeform, integrationorder); }
double parameter::integrate(int physreg, int integrationorder) { return ((expression)*this).integrate(physreg, integrationorder); }

void parameter::write(int physreg, int numfftharms, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, filename, lagrangeorder); }
void parameter::write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder) { return ((expression)*this).write(physreg, numfftharms, meshdeform, filename, lagrangeorder); }

void parameter::write(int physreg, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, filename, lagrangeorder, numtimesteps); }
void parameter::write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps) { return ((expression)*this).write(physreg, meshdeform, filename, lagrangeorder, numtimesteps); }


expression parameter::operator+(void) { return (expression)*this; }
expression parameter::operator-(void) { return -(expression)*this; }

expression parameter::operator+(parameter& inputparameter) { return (expression)*this + inputparameter; }
expression parameter::operator-(parameter& inputparameter) { return (expression)*this - inputparameter; }
expression parameter::operator*(parameter& inputparameter) { return (expression)*this * inputparameter; }
expression parameter::operator/(parameter& inputparameter) { return (expression)*this / inputparameter; }

expression parameter::operator+(double val) { return (expression)*this + val; }
expression parameter::operator-(double val) { return (expression)*this - val; }
expression parameter::operator*(double val) { return (expression)*this * val; }
expression parameter::operator/(double val) { return (expression)*this / val; }


expression operator+(double val, parameter& inputparameter) { return inputparameter+val; }
expression operator-(double val, parameter& inputparameter) { return -inputparameter+val; }
expression operator*(double val, parameter& inputparameter) { return inputparameter*val; }
expression operator/(double val, parameter& inputparameter) { return ( (expression)val )/( (expression)inputparameter ); }





