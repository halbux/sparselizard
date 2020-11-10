#include "rawparameter.h"
#include "disjointregions.h"
#include "disjointregionselector.h"


void rawparameter::synchronize(void)
{
    if (issynchronizing || universe::mymesh->getmeshnumber() == mymeshnumber)
        return;
    issynchronizing = true;    


    // Flush the structure:
    myoperations = std::vector<std::vector<std::shared_ptr<operation>>>(universe::mymesh->getdisjointregions()->count(), std::vector<std::shared_ptr<operation>>(1, NULL));
    maxopnum = -1;
    opnums = std::vector<int>(universe::mymesh->getdisjointregions()->count(),-1);

    // Rebuild the structure:
    for (size_t i = 0; i < mystructuretracker.size(); i++)
        set(mystructuretracker[i].first, mystructuretracker[i].second);
    
    
    mymeshnumber = universe::mymesh->getmeshnumber();
    issynchronizing = false;
}

void rawparameter::errorifundefined(std::vector<int> disjregs)
{
    synchronize();

    for (size_t i = 0; i < disjregs.size(); i++)
    {
        if (myoperations[disjregs[i]].size() == 1 && myoperations[disjregs[i]][0] == NULL)
        {
            std::cout << "Error in 'parameter' object: the parameter has not been defined on the requested region" << std::endl;
            abort();
        }
    }
}

std::vector<int> rawparameter::getopnums(std::vector<int> disjregs)
{
    synchronize();
    
    std::vector<int> output(disjregs.size());
    for (size_t i = 0; i < disjregs.size(); i++)
        output[i] = opnums[disjregs[i]];
    return output;
}

rawparameter::rawparameter(int numrows, int numcols) : myoperations((universe::mymesh->getdisjointregions())->count(), std::vector<std::shared_ptr<operation>>(numrows*numcols, NULL)), opnums((universe::mymesh->getdisjointregions())->count(),-1)
{ 
    mynumrows = numrows; 
    mynumcols = numcols; 
    
    mymeshnumber = universe::mymesh->getmeshnumber();
}

void rawparameter::set(int physreg, expression input)
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

    for (size_t i = 0; i < selecteddisjregs.size(); i++)
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

std::shared_ptr<operation> rawparameter::get(int disjreg, int row, int col)
{
    synchronize();
    
    errorifundefined({disjreg});
    return myoperations[disjreg][row*mynumcols+col];
}

int rawparameter::countrows(void)
{
    return mynumrows;
}

int rawparameter::countcolumns(void)
{
    return mynumcols;
}

std::vector<std::vector<densematrix>> rawparameter::interpolate(int row, int col, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
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
        for (size_t h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1 && out[h].size() == 0)
                out[h] = {densematrix(elemselect.countincurrentorientation(), evaluationcoordinates.size()/3, 0)};
        }
        
        // Insert 'currentinterp' in 'out':
        for (size_t h = 0; h < currentinterp.size(); h++)
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

densematrix rawparameter::multiharmonicinterpolate(int row, int col, int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
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
        for (size_t j = 0; j < selectedelementindexes.size(); j++)
        {
            for (size_t k = 0; k < evaluationcoordinates.size()/3; k++)
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
    
void rawparameter::simplify(int row, int col, int disjreg)
{
    synchronize();
    
    // Make sure the parameter has been defined:
    errorifundefined({disjreg});
    
    myoperations[disjreg][row*mynumcols+col] = myoperations[disjreg][row*mynumcols+col]->simplify({disjreg});
}

void rawparameter::print(void)
{
    synchronize();
    
    std::cout << std::endl;
    std::cout << "Printing parameter of size " << mynumrows << "x" << mynumcols;
    std::cout << std::endl;
    for (size_t i = 0; i < myoperations.size(); i++)
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

