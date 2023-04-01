#include "rawparameter.h"
#include "disjointregions.h"
#include "disjointregionselector.h"


void rawparameter::synchronize(void)
{
    if (issynchronizing || universe::getrawmesh()->getmeshnumber() == mymeshnumber)
        return;
    issynchronizing = true;    


    // Flush the structure:
    myoperations = std::vector<std::vector<std::shared_ptr<operation>>>(universe::getrawmesh()->getdisjointregions()->count(), std::vector<std::shared_ptr<operation>>(mynumrows*mynumcols, NULL));
    maxopnum = -1;
    opnums = std::vector<int>(universe::getrawmesh()->getdisjointregions()->count(),-1);

    // Rebuild the structure:
    for (int i = 0; i < mystructuretracker.size(); i++)
        set(mystructuretracker[i].first, mystructuretracker[i].second);
    
    
    mymeshnumber = universe::getrawmesh()->getmeshnumber();
    issynchronizing = false;
}

void rawparameter::errorifundefined(std::vector<int> disjregs)
{
    synchronize();

    for (int i = 0; i < disjregs.size(); i++)
    {
        if (myoperations[disjregs[i]].size() == 0 || myoperations[disjregs[i]][0] == NULL)
        {
            logs log;
            log.msg() << "Error in 'parameter' object: the parameter has not been defined on the requested region" << std::endl;
            log.error();
        }
    }
}

std::vector<int> rawparameter::getopnums(std::vector<int> disjregs)
{
    synchronize();
    
    std::vector<int> output(disjregs.size());
    for (int i = 0; i < disjregs.size(); i++)
        output[i] = opnums[disjregs[i]];
    return output;
}

rawparameter::rawparameter(int numrows, int numcols) : myoperations(universe::getrawmesh()->getdisjointregions()->count(), std::vector<std::shared_ptr<operation>>(numrows*numcols, NULL)), opnums((universe::getrawmesh()->getdisjointregions())->count(),-1)
{ 
    mynumrows = numrows; 
    mynumcols = numcols; 
    
    mymeshnumber = universe::getrawmesh()->getmeshnumber();
}

void rawparameter::set(int physreg, expression input)
{
    synchronize();
    
    // Keep track of the calls to 'set':
    if (issynchronizing == false)
        mystructuretracker.push_back(std::make_pair(physreg, input));
        
    if (mynumrows != input.countrows() || mynumcols != input.countcolumns())
    {
        logs log;
        log.msg() << "Error in 'parameter' object: trying to set the " << mynumrows << "x" << mynumcols << " sized parameter to a size " << input.countrows() << "x" << input.countcolumns() << std::endl;
        log.error();
    }
    
    maxopnum++;
    
    // Consider ALL disjoint regions in the physical region with (-1):
    std::vector<int> selecteddisjregs = universe::getrawmesh()->getphysicalregions()->get(physreg)->getdisjointregions(-1);

    for (int row = 0; row < mynumrows; row++)
    {
        for (int col = 0; col < mynumcols; col++)
        {
            std::shared_ptr<operation> op = input.getoperationinarray(row, col);
            // Make sure there is no dof or tf in the operation:
            if (op->isdofincluded() || op->istfincluded())
            {
                logs log;
                log.msg() << "Error in 'parameter' object: cannot set an expression containing a dof or a tf" << std::endl;
                log.error();
            }
            // Make sure there is no recursion:
            if (op->isparameterincluded(selecteddisjregs, this))
            {
                logs log;
                log.msg() << "Error in 'parameter' object: cannot set an expression including the parameter itself" << std::endl;
                log.error();
            }
            
            for (int i = 0; i < selecteddisjregs.size(); i++)
            {
                opnums[selecteddisjregs[i]] = maxopnum;
                myoperations[selecteddisjregs[i]][row*mynumcols+col] = op;
            }
        }
    }
}

bool rawparameter::isdefined(int disjreg)
{
    synchronize();
    
    return (myoperations[disjreg].size() > 0 && myoperations[disjreg][0] != NULL);
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

std::vector<std::vector<densemat>> rawparameter::interpolate(int row, int col, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    synchronize();
    
    int numelems = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();
    
    // Make sure the parameter has been defined:
    errorifundefined(alldisjregs);
    
    std::vector<std::vector<densemat>> out = {};
    
    // Group disj. regs. with same operation number (and same element type number).
    disjointregionselector mydisjregselector(alldisjregs, {getopnums(alldisjregs)});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
    
        elemselect.selectdisjointregions(mydisjregs);
        if (elemselect.countinselection() == 0)
            continue;
            
        std::vector<int> selectedelementindexes = elemselect.getelementindexes();
        elementselector myselection = elemselect.extractselection();
        
        auto allstorage = universe::selectsubset(numevalpts, selectedelementindexes);
        // IMPORTANT: Harmonic numbers can be different from one disj. reg. to the other:
        std::vector<std::vector<densemat>> currentinterp = myoperations[mydisjregs[0]][row*mynumcols+col]->interpolate(myselection, evaluationcoordinates, meshdeform);
        universe::restore(allstorage);
        
        // Preallocate the harmonics not yet in 'out':
        if (out.size() < currentinterp.size())
            out.resize(currentinterp.size());
        for (int h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1 && out[h].size() == 0)
                out[h] = {densemat(numelems, numevalpts, 0)};
        }
        
        // Insert 'currentinterp' in 'out':
        for (int h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1)
                out[h][0].insertatrows(selectedelementindexes, currentinterp[h][0]);
        }
    }
    // Unselect the disjoint regions:
    elemselect.selectdisjointregions({});
    
    return out;
}

densemat rawparameter::multiharmonicinterpolate(int row, int col, int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    synchronize();
    
    int numelems = elemselect.countinselection();
    int numevalpts = evaluationcoordinates.size()/3;
    
    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();
    
    // Make sure the parameter has been defined:
    errorifundefined(alldisjregs);
    
    // Preallocate the output matrix:
    densemat out(numtimeevals, numelems * numevalpts);
    
    // Group disj. regs. with same operation number (and same element type number).
    disjointregionselector mydisjregselector(alldisjregs, {getopnums(alldisjregs)});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
        
        elemselect.selectdisjointregions(mydisjregs);
        if (elemselect.countinselection() == 0)
            continue;
            
        std::vector<int> selectedelementindexes = elemselect.getelementindexes();
        elementselector myselection = elemselect.extractselection();
        
        std::vector<int> selectedcolumns(selectedelementindexes.size()*numevalpts);
        for (int j = 0; j < selectedelementindexes.size(); j++)
        {
            for (int k = 0; k < numevalpts; k++)
                selectedcolumns[j*numevalpts+k] = selectedelementindexes[j]*numevalpts+k;
        }
        
        auto allstorage = universe::selectsubset(numevalpts, selectedelementindexes);
        densemat currentinterp = myoperations[mydisjregs[0]][row*mynumcols+col]->multiharmonicinterpolate(numtimeevals, myselection, evaluationcoordinates, meshdeform);
        universe::restore(allstorage);
        
        out.insertatcolumns(selectedcolumns, currentinterp);
    }
    // Unselect the disjoint regions:
    elemselect.selectdisjointregions({});
    
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

