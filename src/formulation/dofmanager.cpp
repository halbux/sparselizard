#include "dofmanager.h"


void dofmanager::synchronize(void)
{
    if (isitmanaged == false || issynchronizing || universe::mymesh->getmeshnumber() == mymeshnumber)
        return;
    issynchronizing = true;    


    // Flush the structure:
    numberofdofs = 0;
    myfields = {};
    int selectedfieldnumberbkp = selectedfieldnumber;
    selectedfieldnumber = -1;
    myfieldorders = {};
    myrawportmap.clear();
    primalondisjreg = {};
    rangebegin = {};
    rangeend = {};

    // Rebuild the structure:
    for (int i = 0; i < myportstructuretracker.size(); i++)
        addtostructure(myportstructuretracker[i]);
    for (int i = 0; i < mystructuretracker.size(); i++)
        addtostructure(mystructuretracker[i].first, mystructuretracker[i].second);
    
    // Select the same field again:
    selectedfieldnumber = selectedfieldnumberbkp;
    
    mymeshnumber = universe::mymesh->getmeshnumber();
    issynchronizing = false;
}

void dofmanager::addtostructure(std::shared_ptr<rawfield> fieldtoadd, std::vector<int> selecteddisjointregions)
{  
    synchronize();
    
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();

    // Find the field index of 'fieldtoadd' (if present):
    int fieldindex = -1;
    for (int i = 0; i < myfields.size(); i++)
    {
        if (myfields[i].get() == fieldtoadd.get())
        {
            fieldindex = i;
            break;
        }
    }
    // Add the field to the structure if not existing.
    // Add an entry for every disjoint region at the same time.
    if (fieldindex == -1)
    {
        fieldindex = myfields.size();
        myfields.push_back(fieldtoadd);
        myfieldorders.push_back(fieldtoadd->getinterpolationorders());
        int numberofdisjointregions = mydisjointregions->count();
        primalondisjreg.push_back(std::vector<std::shared_ptr<rawport>>(numberofdisjointregions, NULL));
        std::vector<std::vector< int >> temp(numberofdisjointregions, std::vector< int >(0));
        rangebegin.push_back(temp);
        rangeend.push_back(temp);
    }
    else
    {
        // Make sure the field order has not changed between two calls:
        if (fieldtoadd->getinterpolationorders() != myfieldorders[fieldindex])
        {
            std::cout << "Error in 'dofmanager' object: ";
            fieldtoadd->print();
            std::cout << " order was changed and does not match dof structure anymore" << std::endl;
            abort();
        }
    }
    
    // Add an entry for every form function - if not already existing:
    for (int i = 0; i < selecteddisjointregions.size(); i++)
    {
        int disjreg = selecteddisjointregions[i];
        
        // Get the element type number in the current disjoint region:
        int elementtypenumber = mydisjointregions->getelementtypenumber(disjreg);
        int elementdimension = mydisjointregions->getelementdimension(disjreg);

        // Get the number of form functions:
        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, fieldtoadd->gettypename());
        int numberofformfunctions = myformfunction->count(fieldtoadd->getinterpolationorder(disjreg), elementdimension, 0);
        
        // Only treat the form functions not yet in the dof structure.
        if (rangebegin[fieldindex][disjreg].size() < numberofformfunctions)
        {
            int numffdefinedbeforeresize = rangebegin[fieldindex][disjreg].size();
            int currentnumberofdofs = mydisjointregions->countelements(disjreg);

            rangebegin[fieldindex][disjreg].resize(numberofformfunctions);
            rangeend[fieldindex][disjreg].resize(numberofformfunctions);
                
            for (int ff = numffdefinedbeforeresize; ff < numberofformfunctions; ff++)
            {
                std::shared_ptr<rawport> pdr = primalondisjreg[fieldindex][disjreg];
                if (pdr != NULL)
                {
                    // Add the primal to the hashmap if not already there:
                    bool isprimalnotthere = (myrawportmap.find(pdr.get()) == myrawportmap.end());
                    if (isprimalnotthere)
                    {
                        myrawportmap[pdr.get()] = numberofdofs;
                        numberofdofs++;
                    }
                    int primaladdress = myrawportmap[pdr.get()];

                    rangebegin[fieldindex][disjreg][ff] = primaladdress;
                    rangeend[fieldindex][disjreg][ff] = primaladdress;
                }
                else
                {
                    rangebegin[fieldindex][disjreg][ff] = numberofdofs;
                    rangeend[fieldindex][disjreg][ff] = numberofdofs + currentnumberofdofs - 1;
                    numberofdofs += currentnumberofdofs;
                }
            }
        }
    }
}

dofmanager::dofmanager(void)
{
    mymeshnumber = universe::mymesh->getmeshnumber();
}

dofmanager::dofmanager(int numdofs)
{
    numberofdofs = numdofs;
    
    isitmanaged = false;
}

void dofmanager::donotsynchronize(void)
{
    issynchronizing = true;
}

void dofmanager::addtostructure(std::shared_ptr<rawport> porttoadd)
{
    synchronize();
    
    // Keep track of the calls to 'addtostructure':
    if (issynchronizing == false)
        myportstructuretracker.push_back(porttoadd);
        
    if (porttoadd->isassociated())
    {
        std::shared_ptr<rawfield> associatedfield = porttoadd->getrawfield();

        disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
        physicalregions* myphysicalregions = universe::mymesh->getphysicalregions();

        // Find the field index of 'associatedfield' (if present):
        int fieldindex = -1;
        for (int i = 0; i < myfields.size(); i++)
        {
            if (myfields[i].get() == associatedfield.get())
            {
                fieldindex = i;
                break;
            }
        }
        // Add the field to the structure if not existing.
        // Add an entry for every disjoint region at the same time.
        if (fieldindex == -1)
        {
            fieldindex = myfields.size();
            myfields.push_back(associatedfield);
            myfieldorders.push_back(associatedfield->getinterpolationorders());
            int numberofdisjointregions = mydisjointregions->count();
            primalondisjreg.push_back(std::vector<std::shared_ptr<rawport>>(numberofdisjointregions, NULL));
            std::vector<std::vector< int >> temp(numberofdisjointregions, std::vector< int >(0));
            rangebegin.push_back(temp);
            rangeend.push_back(temp);
        }
        else
        {
            // Make sure the field order has not changed between two calls:
            if (associatedfield->getinterpolationorders() != myfieldorders[fieldindex])
            {
                std::cout << "Error in 'dofmanager' object: ";
                associatedfield->print();
                std::cout << " order was changed and does not match dof structure anymore" << std::endl;
                abort();
            }
        }

        std::vector<int> disjregs = myphysicalregions->get(porttoadd->getphysicalregion())->getdisjointregions(-1);

        for (int i = 0; i < disjregs.size(); i++)
            primalondisjreg[fieldindex][disjregs[i]] = porttoadd->getprimal();

        // Port to add below to the hashmap:
        porttoadd = porttoadd->getdual();
    }

    // Add the dual to the hashmap (the primal will be added during a regular 'addtostructure' call):
    bool isnotthere = (myrawportmap.find(porttoadd.get()) == myrawportmap.end());
    if (isnotthere)
    {
        myrawportmap[porttoadd.get()] = numberofdofs;
        numberofdofs++;
    }
}

void dofmanager::addtostructure(std::shared_ptr<rawfield> fieldtoadd, int physicalregionnumber)
{
    synchronize();
    
    // Keep track of the calls to 'addtostructure':
    if (issynchronizing == false)
        mystructuretracker.push_back(std::make_pair(fieldtoadd, physicalregionnumber));

    // Get all disjoint regions in the physical region with (-1):
    std::vector<int> disjregs = ((universe::mymesh->getphysicalregions())->get(physicalregionnumber))->getdisjointregions(-1);
    
    addtostructure(fieldtoadd, disjregs);
}

void dofmanager::selectfield(std::shared_ptr<rawfield> selectedfield)
{
    synchronize();
    
    selectedfieldnumber = -1;
    for (int i = 0; i < myfields.size(); i++)
    {
        if (myfields[i].get() == selectedfield.get())
        {
            selectedfieldnumber = i;
            break;
        }
    }
    if (selectedfieldnumber == -1)
    {
        std::cout << "Error in 'dofmanager' object: selected field is not defined in the dof manager" << std::endl;
        abort();
    }
    if (not(issynchronizing) && myfields[selectedfieldnumber]->getinterpolationorders() != myfieldorders[selectedfieldnumber])
    {
        std::cout << "Error in 'dofmanager' object: ";
        myfields[selectedfieldnumber]->print();
        std::cout << " order was changed and does not match dof structure anymore" << std::endl;
        abort();
    }
}

std::vector<int> dofmanager::getdisjointregionsofselectedfield(void)
{
    synchronize();
    
    std::vector<int> output(rangebegin[selectedfieldnumber].size());
    int index = 0;
    for (int i = 0; i < output.size(); i++)
    {
        if (rangebegin[selectedfieldnumber][i].size() > 0)
        {
            output[index] = i;
            index++;
        }
    }
    output.resize(index);
    
    return output;
}

int dofmanager::getrangebegin(int disjreg, int formfunc)
{
    synchronize();
    
    return rangebegin[selectedfieldnumber][disjreg][formfunc];
}

int dofmanager::getrangeend(int disjreg, int formfunc)
{
    synchronize();
    
    return rangeend[selectedfieldnumber][disjreg][formfunc];
}

int dofmanager::getaddress(rawport* prt)
{
    synchronize();
    
    bool isnotthere = (myrawportmap.find(prt) == myrawportmap.end());
    if (isnotthere == false)
        return myrawportmap[prt];
    else
    {
        std::string pn = prt->getname();
        if (pn.size() > 0)
            pn = "'"+pn+"' ";
        std::cout << "Error in 'dofmanager' object: requested port " << pn << "could not be found in the dof structure" << std::endl;
        abort();
    }
}

void dofmanager::getportsinds(std::vector<rawport*>& rps, intdensematrix& inds)
{
    synchronize();
    
    int numports = myrawportmap.size();
    
    rps.resize(numports);
    inds = intdensematrix(numports, 1);
    int* iptr = inds.getvalues();

    int index = 0;
    for (auto it = myrawportmap.begin(); it != myrawportmap.end(); it++)
    {
        rps[index] = it->first;
        iptr[index] = it->second;
        
        index++;
    }
}

std::pair<intdensematrix, intdensematrix> dofmanager::findassociatedports(void)
{
    synchronize();

    intdensematrix primalads(countports(), 1);
    intdensematrix dualads(countports(), 1);

    int* paptr = primalads.getvalues();
    int* daptr = dualads.getvalues();

    int index = 0;
    for (auto it = myrawportmap.begin(); it != myrawportmap.end(); it++)
    {
        rawport* rp = it->first;
        if (rp->isassociated() && rp->isprimal())
        {
            paptr[index] = it->second;
            daptr[index] = myrawportmap[rp->getdual().get()];

            index++;
        }
    }

    primalads = primalads.getresized(index, 1);
    dualads = dualads.getresized(index, 1);

    return std::make_pair(primalads, dualads);
}

bool dofmanager::isdefined(int disjreg, int formfunc)
{
    synchronize();
    
    return (formfunc < rangebegin[selectedfieldnumber][disjreg].size()); 
}

bool dofmanager::isported(int disjreg)
{
    synchronize();
    
    return (primalondisjreg[selectedfieldnumber][disjreg] != NULL); 
}

std::vector<bool> dofmanager::isconstrained(void)
{
    synchronize();
    
    std::vector<bool> output(numberofdofs, false);
    
    std::vector<intdensematrix> allconstrinds = {getdisjregconstrainedindexes(), getgaugedindexes(), getconditionalconstraintdata().first};
    
    for (int i = 0; i < allconstrinds.size(); i++)
    {
        int* cptr = allconstrinds[i].getvalues();
        for (int j = 0; j < allconstrinds[i].count(); j++)
            output[cptr[j]] = true;
    }

    return output;
}

intdensematrix dofmanager::getconstrainedindexes(void)
{
    std::vector<bool> isconstr = isconstrained();

    int numconstr = myalgorithm::counttrue(isconstr);
    intdensematrix out(numconstr, 1);
    int* outptr = out.getvalues();
    
    int index = 0;
    for (int i = 0; i < isconstr.size(); i++)
    {
        if (isconstr[i])
        {
            outptr[index] = i;
            index++;
        }
    }
    
    return out;
}
        
int dofmanager::countdisjregconstraineddofs(void)
{
    synchronize();
    
    int numdisjregconstraineddofs = 0;
    
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        for (int disjreg = 0; disjreg < rangebegin[fieldindex].size(); disjreg++)
        {
            // If the field is constrained on the disjoint region and there is at least one form function:
            if (rangebegin[fieldindex][disjreg].size() > 0 && myfields[fieldindex]->isdisjregconstrained(disjreg))
                numdisjregconstraineddofs += rangebegin[fieldindex][disjreg].size() * (rangeend[fieldindex][disjreg][0] - rangebegin[fieldindex][disjreg][0] + 1);
        }
    }
    return numdisjregconstraineddofs;
}

intdensematrix dofmanager::getdisjregconstrainedindexes(void)
{
    synchronize();
    
    intdensematrix output(1,countdisjregconstraineddofs());
    int* myval = output.getvalues();
    
    int currentindex = 0;
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        for (int disjreg = 0; disjreg < rangebegin[fieldindex].size(); disjreg++)
        {
            // If the field is constrained on the disjoint region.
            if (myfields[fieldindex]->isdisjregconstrained(disjreg))
            {
                for (int ff = 0; ff < rangebegin[fieldindex][disjreg].size(); ff++)
                {
                    int numdofshere = rangeend[fieldindex][disjreg][0] - rangebegin[fieldindex][disjreg][0] + 1;
                    for (int i = 0; i < numdofshere; i++)
                    {
                        myval[currentindex] = rangebegin[fieldindex][disjreg][ff] + i;
                        currentindex++;
                    }
                }
            }
        }
    }
    return output;
}

int dofmanager::countgaugeddofs(void)
{
    synchronize();
    
    int numgaugeddofs = 0;
    
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        for (int disjreg = 0; disjreg < rangebegin[fieldindex].size(); disjreg++)
        {
            if (myfields[fieldindex]->isgauged(disjreg))
            {
                spanningtree* myspantree = myfields[fieldindex]->getspanningtree();

                int elementtype = universe::mymesh->getdisjointregions()->getelementtypenumber(disjreg);
                int fieldorder = myfields[fieldindex]->getinterpolationorder(disjreg);
                std::shared_ptr<hierarchicalformfunction> formfunc = selector::select(elementtype, myfields[fieldindex]->gettypename());
                std::vector<bool> isitgradienttype = formfunc->isgradienttype(fieldorder);
       
                for (int ff = 0; ff < rangebegin[fieldindex][disjreg].size(); ff++)
                {
                    // The lowest order hcurl form function is gauged only on the spanning tree:
                    if (elementtype == 1 && ff == 0)
                        numgaugeddofs += myspantree->countedgesintree(disjreg);
                    else
                    {
                        // All other form functions are gauged on all dofs in the disjoint region:
                        if (isitgradienttype[ff])
                            numgaugeddofs += rangeend[fieldindex][disjreg][0] - rangebegin[fieldindex][disjreg][0] + 1;
                    }
                }
            }
        }
    }
    return numgaugeddofs;
}

intdensematrix dofmanager::getgaugedindexes(void)
{
    synchronize();
    
    intdensematrix output(1,countgaugeddofs());
    int* myval = output.getvalues();
    
    int currentindex = 0;
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        for (int disjreg = 0; disjreg < rangebegin[fieldindex].size(); disjreg++)
        {
            if (myfields[fieldindex]->isgauged(disjreg))
            {
                spanningtree* myspantree = myfields[fieldindex]->getspanningtree();

                int elementtype = universe::mymesh->getdisjointregions()->getelementtypenumber(disjreg);
                int fieldorder = myfields[fieldindex]->getinterpolationorder(disjreg);
                std::shared_ptr<hierarchicalformfunction> formfunc = selector::select(elementtype, myfields[fieldindex]->gettypename());
                std::vector<bool> isitgradienttype = formfunc->isgradienttype(fieldorder);
       
                for (int ff = 0; ff < rangebegin[fieldindex][disjreg].size(); ff++)
                {
                    // The lowest order hcurl form function is gauged only on the spanning tree.
                    // All other form functions are gauged on all dofs in the disjoint region.
                    int numdofshere = rangeend[fieldindex][disjreg][0] - rangebegin[fieldindex][disjreg][0] + 1;
                    for (int i = 0; i < numdofshere; i++)
                    {
                        if ((elementtype != 1 || ff != 0) && isitgradienttype[ff] || elementtype == 1 && ff == 0 && myspantree->isintree(i, disjreg))
                        {
                            myval[currentindex] = rangebegin[fieldindex][disjreg][ff] + i;
                            currentindex++;
                        }
                    }
                }
            }
        }
    }
    return output;
}

std::pair<intdensematrix, densematrix> dofmanager::getconditionalconstraintdata(void)
{
    synchronize();
    
    // This will have an entry for every field and every disjoint node region that is conditionally constrained:
    std::vector<intdensematrix> indexmat = {};
    std::vector<densematrix> condvalvec = {};
    std::vector<densematrix> constrvalvec = {};
    
    
    for (int fieldindex = 0; fieldindex < rangebegin.size(); fieldindex++)
    {
        // First get the list of disjoint NODE regions on which the rawfield is conditionally constrained:
        std::vector<bool> isdisjregactive(rangebegin[fieldindex].size(), false);
        for (int disjreg = 0; disjreg < isdisjregactive.size(); disjreg++)
        {
            // Only the nodes are constrained:
            if (rangebegin[fieldindex][disjreg].size() == 0 || universe::mymesh->getdisjointregions()->getelementtypenumber(disjreg) != 0)
                continue;
                
            if (myfields[fieldindex]->isconditionallyconstrained(disjreg))
                isdisjregactive[disjreg] = true;
        }
            
        // Get the condition and value expressions for the conditional constraint:
        std::vector<std::vector<expression>> condconstrexpr = myfields[fieldindex]->getconditionalconstraints();
    
        // Loop on all disjoint regions:
        for (int disjreg = 0; disjreg < isdisjregactive.size(); disjreg++)
        {
            if (isdisjregactive[disjreg] == false)
                continue;
            
            std::shared_ptr<operation> condop = condconstrexpr[disjreg][0].getoperationinarray(0,0), constrop = condconstrexpr[disjreg][1].getoperationinarray(0,0);
            
            // Combine all disjoint regions that share the same condop and constrop operations:
            std::vector<int> curdisjregs = {};
            for (int i = disjreg; i < isdisjregactive.size(); i++)
            {
                if (isdisjregactive[i] && condop == condconstrexpr[i][0].getoperationinarray(0,0) && constrop == condconstrexpr[i][1].getoperationinarray(0,0))
                {
                    curdisjregs.push_back(i);
                    isdisjregactive[i] = false;
                }
            }
                
            // Compute the conditional and constraint expressions.
            // For nodes the reference coordinates are all zero.
            std::vector<double> evaluationcoordinates = {0,0,0};
            elementselector myelemselect(curdisjregs, condop->isvalueorientationdependent(curdisjregs) || constrop->isvalueorientationdependent(curdisjregs));

            std::vector<std::vector<densematrix>> condval = condop->interpolate(myelemselect, evaluationcoordinates, NULL);
            // Skip if no conditional constraint is active:
            if (condval[1][0].max() < 0)
                continue;
            std::vector<std::vector<densematrix>> constrval = constrop->interpolate(myelemselect, evaluationcoordinates, NULL);

            // Create a matrix with all indices:
            intdensematrix curindexmat(condval[1][0].countrows(), condval[1][0].countcolumns(),0);
            int* curindexmatptr = curindexmat.getvalues();
            int index = 0;
            for (int i = 0; i < curdisjregs.size(); i++)
            {
                // There is only a single shape function per node!
                for (int ind = rangebegin[fieldindex][curdisjregs[i]][0]; ind <= rangeend[fieldindex][curdisjregs[i]][0]; ind++)
                {
                    curindexmatptr[index] = ind;
                    index++;
                }
            }

            // Append to the vectors:
            condvalvec.push_back(condval[1][0]);
            constrvalvec.push_back(constrval[1][0]);
            indexmat.push_back(curindexmat);
        }
    }

    // Count the number of active conditional constraints:
    int numactive = 0;
    for (int i = 0; i < condvalvec.size(); i++)
    {
        double* condvalptr = condvalvec[i].getvalues();
        for (int j = 0; j < condvalvec[i].count(); j++)
        {
            if (condvalptr[j] >= 0)
                numactive++;
        }
    }
    
    // Combine all active conditional constraints together:
    intdensematrix condconstrindices(numactive,1);
    densematrix condconstrval(numactive,1);
    
    int* indptr = condconstrindices.getvalues();
    double* valptr = condconstrval.getvalues();
    
    int index = 0;
    for (int i = 0; i < condvalvec.size(); i++)
    {
        double* condvalptr = condvalvec[i].getvalues();
        double* constrvalptr = constrvalvec[i].getvalues();
        int* indmatptr = indexmat[i].getvalues();
        for (int j = 0; j < condvalvec[i].count(); j++)
        {
            if (condvalptr[j] >= 0)
            {
                indptr[index] = indmatptr[j];
                valptr[index] = constrvalptr[j];
                
                index++;
            }
        }
    }
    
    return std::make_pair(condconstrindices, condconstrval);
}

std::shared_ptr<rawfield> dofmanager::getselectedfield(void)
{
    synchronize();

    if (selectedfieldnumber >= 0)
        return myfields[selectedfieldnumber];
    else
        return NULL;
}

std::vector<std::shared_ptr<rawfield>> dofmanager::getfields(void)
{
    synchronize();
    
    return myfields;
}

void dofmanager::replaceselectedfield(std::shared_ptr<rawfield> rf)
{
    synchronize();
    
    myfields[selectedfieldnumber] = rf;
}

std::vector<int> dofmanager::getselectedfieldorders(void)
{
    synchronize();
    
    return myfieldorders[selectedfieldnumber];
}

int dofmanager::countports(void)
{
    synchronize();

    return myrawportmap.size();
}

int dofmanager::countassociatedprimalports(void)
{
    synchronize();

    int cnt = 0;
    for (auto it = myrawportmap.begin(); it != myrawportmap.end(); it++)
    {
        rawport* rp = it->first;
        if (rp->isassociated() && rp->isprimal())
            cnt++;
    }

    return cnt;
}

int dofmanager::countdofs(void)
{
    synchronize();
    
    return numberofdofs;
}

long long int dofmanager::allcountdofs(void)
{
    synchronize();
    
    if (slmpi::count() == 1)
        return countdofs();
    
    std::shared_ptr<dtracker> mydtracker = universe::mymesh->getdtracker();
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
    
    mydtracker->errorundefined();
    
    long long int out = 0;
    std::vector<bool> isowndr = mydtracker->isdisjointregionowned();
    
    for (int d = 0; d < isowndr.size(); d++)
    {
        if (isowndr[d])
        {
            int ne = mydisjointregions->countelements(d);
            for (int i = 0; i < myfields.size(); i++)
            {
                int nff = rangebegin[i][d].size();
                out += ne*nff;
            }
        }
    }
    
    slmpi::sum(1, &out);
    
    return out;
}

int dofmanager::countformfunctions(int disjointregion)
{
    synchronize();
    
    return rangebegin[selectedfieldnumber][disjointregion].size();
}

std::vector<std::vector<intdensematrix>> dofmanager::discovernewconstraints(std::vector<int> neighbours, std::vector<intdensematrix> senddofinds, std::vector<intdensematrix> recvdofinds)
{
    // New constraints can only appear at the outer overlap/no-overlap interfaces.

    synchronize();
    
    int numneighbours = neighbours.size();
    
    std::vector<intdensematrix> sendnewconstrainedinds(numneighbours), recvnewconstrainedinds(numneighbours), sendunconstrainedinds(numneighbours), recvunconstrainedinds(numneighbours);
    
    // Get all types of constraints:
    std::vector<bool> isdofconstrained = isconstrained(); 
    
    // Exchange the constrained indexes:
    std::vector<std::vector<int>> isconstrainedforneighbours(numneighbours), isconstrainedfromneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        std::vector<bool> ic;
        myalgorithm::select(isdofconstrained, senddofinds[n], ic);

        myalgorithm::pack(ic, isconstrainedforneighbours[n]);
        isconstrainedfromneighbours[n].resize(myalgorithm::getpackedsize(recvdofinds[n].count()));
    }
    
    slmpi::exchange(neighbours, isconstrainedforneighbours, isconstrainedfromneighbours);
    
    std::vector<std::vector<int>> isnewlyconstrainedforneighbours(numneighbours), isnewlyconstrainedfromneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        int* recvindsptr = recvdofinds[n].getvalues();
                
        std::vector<bool> isconstrainedinneighbour;
        myalgorithm::unpack(recvdofinds[n].count(), isconstrainedfromneighbours[n], isconstrainedinneighbour);
        
        std::vector<bool> isnewlyconstrained(recvdofinds[n].count(), false);
        
        for (int i = 0; i < isconstrainedinneighbour.size(); i++)
        {
            int recvind = recvindsptr[i];
            if (isconstrainedinneighbour[i] == true && isdofconstrained[recvind] == false)
            {
                isdofconstrained[recvind] = true;
                isnewlyconstrained[i] = true;
            }
        }
        recvnewconstrainedinds[n] = recvdofinds[n].select(isnewlyconstrained, true);
        
        myalgorithm::pack(isnewlyconstrained, isnewlyconstrainedforneighbours[n]);
        isnewlyconstrainedfromneighbours[n].resize(myalgorithm::getpackedsize(senddofinds[n].count()));
    }
    slmpi::exchange(neighbours, isnewlyconstrainedforneighbours, isnewlyconstrainedfromneighbours);
    
    // List the indexes on this rank of the new constraints in each neighbour:
    for (int n = 0; n < numneighbours; n++)
    {
        std::vector<bool> isnewlyconstrainedinneighbour;
        myalgorithm::unpack(senddofinds[n].count(), isnewlyconstrainedfromneighbours[n], isnewlyconstrainedinneighbour);
        sendnewconstrainedinds[n] = senddofinds[n].select(isnewlyconstrainedinneighbour, true);
    }
    
    // List the unconstrained indexes to exchange with each neighbour:
    for (int n = 0; n < numneighbours; n++)
    {
        std::vector<bool> si, ri;
        myalgorithm::select(isdofconstrained, senddofinds[n], si);
        sendunconstrainedinds[n] = senddofinds[n].select(si, false);
        myalgorithm::select(isdofconstrained, recvdofinds[n], ri);
        recvunconstrainedinds[n] = recvdofinds[n].select(ri, false);
    }
    
    return {sendnewconstrainedinds, recvnewconstrainedinds, sendunconstrainedinds, recvunconstrainedinds};
}
        
void dofmanager::print(void)
{
    synchronize();
    
    std::cout << "Showing the content of the dof manager (" << numberofdofs << " dofs in total):" << std::endl << std::endl;
    
    for (int i = 0; i < myfields.size(); i++)
    {
        for (int disjreg = 0; disjreg < rangebegin[i].size(); disjreg++)
        {
            // Only print not-empty entries:
            if (rangebegin[i][disjreg].size() != 0)
            {
                std::cout << "Field ";
                myfields[i]->print();
                std::cout << " @ disj. reg. " << disjreg << " is in range " << rangebegin[i][disjreg][0] << "-" << rangeend[i][disjreg][rangeend[i][disjreg].size()-1] << std::endl;
            }
        }
    }
    std::cout << std::endl;
}


intdensematrix dofmanager::getaddresses(std::shared_ptr<rawfield> inputfield, int fieldinterpolationorder, int elementtypenumber, std::vector<int> &elementlist, int fieldphysreg)
{
    synchronize();
    
    elements* myelements = universe::mymesh->getelements();
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();

    // Create a vector whose ith index is true if the field is defined 
    // on the disjoint region number i and false otherwise.
    std::vector<int> fielddisjregs = ((universe::mymesh->getphysicalregions())->get(fieldphysreg))->getdisjointregions(-1); // Get all disj regs with (-1)
    std::vector<bool> isfielddefinedondisjointregion(mydisjointregions->count(),false);
    for (int i = 0; i < fielddisjregs.size(); i++)
        isfielddefinedondisjointregion[fielddisjregs[i]] = true;
                    
    element myelement(elementtypenumber);
    
    // Create a form function iterator to iterate through all form functions.
    hierarchicalformfunctioniterator myiterator(inputfield->gettypename(), elementtypenumber, fieldinterpolationorder);

    int numrows = myiterator.count();
    int numcols = elementlist.size();
    
    intdensematrix output(numrows, numcols);
    int* adresses = output.getvalues();

    selectfield(inputfield);
    
    for (int ff = 0; ff < myiterator.count(); ff++)
    {
        int associatedelementtype = myiterator.getassociatedelementtype();
        int associatedelementdimension = myiterator.getdimension();
        int formfunctionindex = myiterator.getformfunctionindexinnodeedgefacevolume();
        int num = myiterator.getnodeedgefacevolumeindex();
        // For quad subelements in prisms and pyramids:
        if ((elementtypenumber == 6 || elementtypenumber == 7) && associatedelementtype == 3)
            num -= myelement.counttriangularfaces();
        
        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, inputfield->gettypename());
        int currentnumberofformfunctions, previousdisjreg;       

        for (int i = 0; i < elementlist.size(); i++)
        {
            int elem = elementlist[i];
            // Get the subelement to which the current form function is associated:
            int currentsubelem = myelements->getsubelement(associatedelementtype, elementtypenumber, elem, num);
            // Also get its disjoint region number:
            int currentdisjointregion = myelements->getdisjointregion(associatedelementtype, currentsubelem);
            
            // The current subelement might require less form functions:
            if (i == 0 || currentdisjointregion != previousdisjreg)
                currentnumberofformfunctions = myformfunction->count(inputfield->getinterpolationorder(currentdisjointregion), associatedelementdimension, myiterator.getnodeedgefacevolumeindex());

            previousdisjreg = currentdisjointregion;
            
            // If not in a disjoint region on which the field is defined set -1 address.
            if (isfielddefinedondisjointregion[currentdisjointregion] && formfunctionindex < currentnumberofformfunctions)
            {
                // Use it to get the subelem index in the disjoint region:
                currentsubelem -= mydisjointregions->getrangebegin(currentdisjointregion);
                
                if (primalondisjreg[selectedfieldnumber][currentdisjointregion] == NULL)
                    adresses[ff*numcols+i] = rangebegin[selectedfieldnumber][currentdisjointregion][formfunctionindex] + currentsubelem;
                else
                    adresses[ff*numcols+i] = rangebegin[selectedfieldnumber][currentdisjointregion][formfunctionindex] + 0;
            }
            else
                adresses[ff*numcols+i] = -1;
        }
        myiterator.next();
    }

    return output;
}

