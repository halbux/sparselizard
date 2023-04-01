#include "rawfield.h"


void rawfield::synchronize(std::vector<int> physregsfororder, std::vector<int> disjregsfororder)
{
    // The coordinate fields can be used even before the mesh is loaded (cannot call the ptracker).
    if (mytypename == "x" || mytypename == "y" || mytypename == "z" || mytypename == "")
        return;
        
    // Synchronize all subfields and harmonics:
    if (mysubfields.size() != 0 || myharmonics.size() != 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->synchronize();
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                myharmonics[h][0]->synchronize();
        }
        return;
    }
        
    if (issynchronizing || not(issynchronizingallowed) || myptracker == universe::getrawmesh()->getptracker())
        return;
    issynchronizing = true; 
    
    // Get a copy of this field before syncing:
    std::shared_ptr<rawfield> originalthis(new rawfield);
    *originalthis = *this;
    
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    
    // Create a new coef manager:
    mycoefmanager = std::shared_ptr<coefmanager>(new coefmanager(mytypename, drs));
    
    // Flush the containers:
    interpolationorder = std::vector<int>(drs->count(), -1);
    mydisjregconstraints = std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>>(drs->count(), NULL);
    myconditionalconstraints = std::vector<std::vector<expression>>(drs->count(), std::vector<expression>(0));
    isitgauged = std::vector<bool>(drs->count(), false);
    isitported = std::vector<bool>(drs->count(), false);
    
    // Rebuild the containers:
    if (disjregsfororder.size() == 0)
    {
        if (physregsfororder.size() == 0)
        {
            for (int i = 0; i < myordertracker.size(); i++)
                setorder(myordertracker[i].first, myordertracker[i].second, false);
        }
        else
        {
            // For p-adaptive fields:
            for (int i = physregsfororder.size()/2-1; i >= 0; i--)
                setorder(physregsfororder[2*i+0], physregsfororder[2*i+1], false);
        }
    }
    else
    {
        interpolationorder = disjregsfororder;
        for (int i = 0; i < interpolationorder.size(); i++)
            mycoefmanager->fitinterpolationorder(i, interpolationorder[i]); 
    }
    for (int i = 0; i < mydisjregconstrainttracker.size(); i++)
    {
        expression* meshdef = NULL;
        if (std::get<2>(mydisjregconstrainttracker[i]).size() > 0)
            meshdef = &(std::get<2>(mydisjregconstrainttracker[i])[0]);
        setdisjregconstraint(std::get<0>(mydisjregconstrainttracker[i]), std::get<1>(mydisjregconstrainttracker[i]), meshdef, std::get<3>(mydisjregconstrainttracker[i]), std::get<4>(mydisjregconstrainttracker[i]));
    }
    for (int i = 0; i < myconditionalconstrainttracker.size(); i++)
        setconditionalconstraint(std::get<0>(myconditionalconstrainttracker[i]), std::get<1>(myconditionalconstrainttracker[i]), std::get<2>(myconditionalconstrainttracker[i]));
    for (int i = 0; i < mygaugetracker.size(); i++)
        setgauge(mygaugetracker[i]);
    // This forces the lowest order on all ports (overwrites also p-adaptivity):
    for (int i = 0; i < myporttracker.size(); i++)
        setport(std::get<0>(myporttracker[i]), std::get<1>(myporttracker[i]), std::get<2>(myporttracker[i]));
    
    // Update the coef manager with the new nodal/edge/face/volume shape function coefficients:
    if (isvaluesynchronizingallowed)
        updateshapefunctions(sl::athp(field(originalthis), myrawmesh, myptracker), NULL, {drs->getindim(0), drs->getindim(1), drs->getindim(2), drs->getindim(3)}, myupdateaccuracy);
    
        
    // Update the mesh tracker to the current one:
    myptracker = universe::getrawmesh()->getptracker();
    myrawmesh = universe::getrawmesh();
    issynchronizing = false;
}

void rawfield::updateshapefunctions(expression updateexpr, expression* meshdeform, std::vector<std::vector<int>> drsindims, int updateaccuracy, bool withtiming)
{
    // Only do a projection. No constraint calculation.
    std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>> mydisjregconstraintsbkp = mydisjregconstraints;
    std::vector<std::vector<expression>> myconditionalconstraintsbkp = myconditionalconstraints;
    std::vector<bool> isitgaugedbkp = isitgauged;
    std::vector<bool> isitportedbkp = isitported;
    
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    
    mydisjregconstraints = std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>>(drs->count(), NULL);
    myconditionalconstraints = std::vector<std::vector<expression>>(drs->count(), std::vector<expression>(0));
    isitgauged = std::vector<bool>(drs->count(), false);
    isitported = std::vector<bool>(drs->count(), false);


    wallclock clkn;
    updatenodalshapefunctions(updateexpr, meshdeform, drsindims);
    if (withtiming)
        clkn.print("Time to update the nodal shape functions:");
    wallclock clke;
    updateothershapefunctions(1, updateexpr, meshdeform, drsindims, updateaccuracy);
    if (withtiming)
        clke.print("Time to update the edge shape functions:");
    wallclock clkf;
    updateothershapefunctions(2, updateexpr, meshdeform, drsindims, updateaccuracy);
    if (withtiming)
        clkf.print("Time to update the face shape functions:");
    wallclock clkv;
    updateothershapefunctions(3, updateexpr, meshdeform, drsindims, updateaccuracy);
    if (withtiming)
    {
        clkv.print("Time to update the volume shape functions:");
        clkn.print("Total time:");
    }
    
    
    // Restore:
    mydisjregconstraints = mydisjregconstraintsbkp;
    myconditionalconstraints = myconditionalconstraintsbkp;
    isitgauged = isitgaugedbkp;
    isitported = isitportedbkp;
}

void rawfield::updatenodalshapefunctions(expression updateexpr, expression* meshdeform, std::vector<std::vector<int>> drsindims)
{
    field thisfield = field(shared_from_this());
    
    if (gettypename() != "h1" || drsindims[0].size() == 0)
        return;
        
    physicalregions* prs = universe::getrawmesh()->getphysicalregions();
        
    // Create temporary physical regions:
    int physreg = prs->createfromdisjointregionlist(0, drsindims[0]);
    
    formulation evalatnodes;
    
    integration myterm;
    if (meshdeform == NULL)
        myterm = integration(physreg, -updateexpr * sl::tf(thisfield));
    else
        myterm = integration(physreg, *meshdeform, -updateexpr * sl::tf(thisfield));
    myterm.isbarycentereval = true;    
    
    evalatnodes += myterm;

    evalatnodes.generaterhs();
    vec vals = evalatnodes.rhs(false, false);
    
    setdata(physreg, vals|thisfield);
    prs->remove({physreg}, false);
}

void rawfield::updateothershapefunctions(int dim, expression updateexpr, expression* meshdeform, std::vector<std::vector<int>> drsindims, int updateaccuracy) // dim can be 1, 2 or 3
{
    if (drsindims[dim].size() == 0)
        return;

    std::string tn = gettypename();
    field thisfield = field(shared_from_this());
    
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    physicalregions* prs = universe::getrawmesh()->getphysicalregions();
    
    // Create temporary physical regions:
    int physreg = prs->createfromdisjointregionlist(dim, drsindims[dim]);
    int dirichletphysreg;
    if (dim == 1)
        dirichletphysreg = prs->createfromdisjointregionlist(0, drsindims[0]);
    if (dim == 2)
        dirichletphysreg = prs->createfromdisjointregionlist(1, gentools::concatenate({drsindims[0],drsindims[1]}));
    if (dim == 3)
        dirichletphysreg = prs->createfromdisjointregionlist(2, gentools::concatenate({drsindims[0],drsindims[1],drsindims[2]}));
    
    // The Dirichlet constraints (if any) are added to the rhs of the projection.
    //
    // | A   B |  | x |   | v |
    // |       |  |   | = |   |
    // | 0   1 |  | y |   | w |

    // Create the formulation to get block A and v:
    formulation blockAv;
    // A and v blocks of the projection:
    if (meshdeform == NULL)
        blockAv += sl::integral(physreg, sl::dof(thisfield) * sl::tf(thisfield) - updateexpr * sl::tf(thisfield), updateaccuracy);
    else
        blockAv += sl::integral(physreg, *meshdeform, sl::dof(thisfield) * sl::tf(thisfield) - updateexpr * sl::tf(thisfield), updateaccuracy);
  
    std::shared_ptr<dofmanager> dm = blockAv.getdofmanager();
    dm->selectfield(shared_from_this());
    
    // Get the block diagonal info:
    std::vector<int> alldrsindim = drsindims[dim];
    // Count the number of non-empty diagonal blocks:
    int numblocks = 0, preallocsize = 0;
    for (int d = 0; d < alldrsindim.size(); d++)
    {
        int curdr = alldrsindim[d];
        int numelemsindr = drs->countelements(curdr);
        int numffindr = dm->countformfunctions(curdr);
        if (numffindr == 0)
            continue;
        numblocks += numelemsindr;
        preallocsize += numelemsindr * numffindr*numffindr;
    }
    
    if (preallocsize > 0)
    {
        indexmat blocksizes(numblocks,1);
        int* bsvals = blocksizes.getvalues();

        // Vector to reorder the mat and vec to bring together the parts of the diagonal blocks:
        indexmat renumtodiagblocks(dm->countdofs(), 1);
        int* renumptr = renumtodiagblocks.getvalues();
        int index = 0;
        for (int d = 0; d < alldrsindim.size(); d++)
        {
            int curdr = alldrsindim[d];
            int numelemsindr = drs->countelements(curdr);
            int numffindr = dm->countformfunctions(curdr);
            if (numffindr == 0)
                continue;
            int rb = dm->getrangebegin(curdr, 0);
            
            for (int i = 0; i < numelemsindr; i++)
            {
                bsvals[index+i] = numffindr;
                for (int ff = 0; ff < numffindr; ff++)
                    renumptr[rb+i*numffindr+ff] = rb+numelemsindr*ff+i;
            }
            index += numelemsindr;
        }
    
        // Get blocks A and v:
        blockAv.generate();
        mat A = blockAv.A();
        vec v = blockAv.rhs(false, false);
        
        if (dim > 1 || tn == "h1")
        {
            // Create the formulation to get block B and w:
            formulation blockB;
            // B block of the projection:
            if (meshdeform == NULL)
                blockB += sl::integral(physreg, sl::dof(thisfield, dirichletphysreg) * sl::tf(thisfield), updateaccuracy);
            else
                blockB += sl::integral(physreg, *meshdeform, sl::dof(thisfield, dirichletphysreg) * sl::tf(thisfield), updateaccuracy);
            // Get block B:
            blockB.generatestiffnessmatrix();
            mat B = blockB.A();
            vec w = vec(blockB);
            transferdata(dirichletphysreg, w|thisfield, "set");
            
            // System to solve is Ax = v - By:
            vec By = -B*w;
            setdata(physreg, By|thisfield);
            transferdata(physreg, v|thisfield, "add");
        }
        
        // Permute A:
        Mat permutedmat;
        IS permutis;
        ISCreateGeneral(PETSC_COMM_SELF, renumtodiagblocks.count(), renumtodiagblocks.getvalues(), PETSC_USE_POINTER, &permutis);
        ISSetPermutation(permutis);
        MatPermute(A.getapetsc(), permutis, permutis, &permutedmat);
        
        densemat blockvals(preallocsize, 1);
        MatInvertVariableBlockDiagonal(permutedmat, blocksizes.count(), bsvals, blockvals.getvalues());
        MatDestroy(&permutedmat);
        
        // Solve block-diagonal system:
        v.permute(renumtodiagblocks);
        indexmat alladds(v.size(),1, 0,1);
        densemat vmat = v.getvalues(alladds);
        densemat prod = blockvals.blockdiagonaltimesvector(blocksizes, vmat);

        v.setvalues(alladds, prod);
        v.permute(renumtodiagblocks, true);

        setdata(physreg, v|thisfield);
    }
    prs->remove({dirichletphysreg}, false);
    prs->remove({physreg}, false);
}

void rawfield::allowsynchronizing(bool allowit)
{
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->allowsynchronizing(allowit);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->allowsynchronizing(allowit);
    }
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
        issynchronizingallowed = allowit;
}

void rawfield::allowvaluesynchronizing(bool allowit)
{
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->allowvaluesynchronizing(allowit);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->allowvaluesynchronizing(allowit);
    }
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
        isvaluesynchronizingallowed = allowit;
}

void rawfield::setupdateaccuracy(int extraintegrationorder)
{
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setupdateaccuracy(extraintegrationorder);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setupdateaccuracy(extraintegrationorder);
    }
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
        myupdateaccuracy = extraintegrationorder;
}

rawfield::rawfield(std::string fieldtypename, const std::vector<int> harmonicnumbers, bool ismultiharm)
{
    amimultiharmonic = ismultiharm;
    
    // Treat the coordinate fields:
    if (fieldtypename == "x" || fieldtypename == "y" || fieldtypename == "z")
    {            
        myname = fieldtypename;
        mytypename = fieldtypename;
        // A coordinate field can only be constant (i.e. on harmonic 1). 
        if (harmonicnumbers.size() != 1 || harmonicnumbers[0] != 1)
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: rawfield type " << fieldtypename << " cannot be multiharmonic (must be constant)" << std::endl;
            log.error();
        }
        return;
    }

    // Make sure the mesh has been loaded before defining non-coordinate fields:
    if (universe::myrawmesh == NULL)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: first load mesh before defining a field that is not the x, y or z coordinate" << std::endl;
        log.error();
    }

    // If the field type name ends with xy (or xyz) there are 2 (or 3) dof components.
    // 'xy' or 'xyz' can only be used on scalar form function type (e.g. not on hcurl).
    mytypename = fieldtypename;

    int numberofsubfields = 1;
    if (mytypename.size() > 2 && mytypename.compare(mytypename.size()-2,2,"xy") == 0)
        numberofsubfields = 2;
    if (mytypename.size() > 3 && mytypename.compare(mytypename.size()-3,3,"xyz") == 0)
        numberofsubfields = 3;
    // Erase any trailing xy or xyz:
    if (numberofsubfields > 1)
        mytypename.erase(mytypename.size()-numberofsubfields,numberofsubfields);  
        
    // Translate to the actual type name:
    if (mytypename == "one")
        mytypename = "one"+std::to_string(universe::getrawmesh()->getmeshdimension());
    if (mytypename == "h1d")
        mytypename = "h1d"+std::to_string(universe::getrawmesh()->getmeshdimension());

    // Make sure the subfield type is scalar:
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(0, mytypename);
    if (numberofsubfields > 1 && myformfunction->countcomponents() != 1)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: " << fieldtypename << " is not a valid field type. Only scalar form functions can have a trailing 'xy' or 'xyz'" << std::endl;
        log.error();
    }

    // For fields with subfields:
    if (numberofsubfields > 1)
    {
        for (int i = 0; i < numberofsubfields; i++)
            mysubfields.push_back({ std::shared_ptr<rawfield>(new rawfield(mytypename, harmonicnumbers, ismultiharm)) });
    }
    else
    {
        if (harmonicnumbers.size() != 0)
        {
            myharmonics.resize(*max_element(harmonicnumbers.begin(), harmonicnumbers.end())+1);
            for (int h = 0; h < harmonicnumbers.size(); h++)
                myharmonics[harmonicnumbers[h]] = { std::shared_ptr<rawfield>(new rawfield(mytypename, {}, ismultiharm)) };
        }    
        else
        {
            myrawmesh = universe::getrawmesh();
        
            myptracker = universe::getrawmesh()->getptracker();
            
            disjointregions* drs = myrawmesh->getdisjointregions();
        
            // Set a -1 undefined interpolation order by default:
            interpolationorder = std::vector<int>(drs->count(), -1);
            // Set all unconstrained by default:
            mydisjregconstraints = std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>>(drs->count(), NULL);
            
            myconditionalconstraints = std::vector<std::vector<expression>>(drs->count(), std::vector<expression>(0));

            isitgauged = std::vector<bool>(drs->count(), false);
            
            isitported = std::vector<bool>(drs->count(), false);

            mycoefmanager = std::shared_ptr<coefmanager>(new coefmanager(mytypename, drs));
        }
    }
    return;
}

rawfield::rawfield(void) {}

rawfield::rawfield(dofmanager* dm, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt)
{
    int numdrs = pt->getdisjointregions()->count();
    
    std::shared_ptr<rawfield> selectedrf = dm->getselectedfield();

    amimultiharmonic = selectedrf->amimultiharmonic;
    mytypename = selectedrf->gettypename();
    mycoefmanager = std::shared_ptr<coefmanager>(new coefmanager(mytypename, pt->getdisjointregions()));
    
    interpolationorder = dm->getselectedfieldorders();
    
    mydisjregconstraints = std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>>(numdrs, NULL);
    myconditionalconstraints = std::vector<std::vector<expression>>(numdrs, std::vector<expression>(0));
    isitgauged = std::vector<bool>(numdrs, false);
    isitported = std::vector<bool>(numdrs, false);
    
    for (int i = 0; i < numdrs; i++)
        mycoefmanager->fitinterpolationorder(i, interpolationorder[i]); 

    issynchronizingallowed = selectedrf->issynchronizingallowed;
    isvaluesynchronizingallowed = selectedrf->isvaluesynchronizingallowed;
    myupdateaccuracy = selectedrf->myupdateaccuracy;
    
    myptracker = pt;
    myrawmesh = rm;
}

rawfield::~rawfield(void)
{
    if (universe::myrawmesh != NULL && mysubfields.size() == 0 && myharmonics.size() == 0)
        universe::getrawmesh()->remove(this);
}

int rawfield::countcomponents(void)
{
    // Do not sync this.

    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
        return 1;
        
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(0, mytypename);
    return std::max(myformfunction->countcomponents(), (int) mysubfields.size());
}

int rawfield::countformfunctioncomponents(void)
{
    // Do not sync this.
    
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
        return 1;
    
    // Get the number of form function components:
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(0, mytypename);
    return myformfunction->countcomponents();
}

std::vector<int> rawfield::getharmonics(void)
{
    synchronize();
    
    if (mysubfields.size() == 0)
    {
        std::vector<int> harms = {};

        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                harms.push_back(h);
        }
        if (myharmonics.size() == 0)
            harms = {1};

        return harms;
    }
    else
        return mysubfields[0][0]->getharmonics();
}

int rawfield::getfirstharmonic(void)
{
    // Do not sync this.
    
    if (mysubfields.size() == 0)
    {
        if (myharmonics.size() == 0)
            return 1;
        else
        {
            for (int h = 0; h < myharmonics.size(); h++)
            {
                if (myharmonics[h].size() == 1)
                    return h;
            }
        }
    }
    else
        return mysubfields[0][0]->getfirstharmonic();
        
    throw std::runtime_error(""); // fix return warning
}

bool rawfield::isharmonicincluded(int harmonic)
{
    synchronize();
    
    if (mysubfields.size() == 0)
        return (harmonic == 1 && myharmonics.size() == 0) || (harmonic > 0 && harmonic < myharmonics.size() && myharmonics[harmonic].size() > 0);
    else
        return mysubfields[0][0]->isharmonicincluded(harmonic);
}

void rawfield::printharmonics(void)
{
    synchronize();
    
    if (amimultiharmonic == false)
    {
        std::cout << "Field is not multiharmonic" << std::endl;
        return;
    }

    if (mysubfields.size() > 0)
        mysubfields[0][0]->printharmonics();
    else
    {
        if (myharmonics.size() == 0)
        {
            std::cout << " +vc0*cos(0*pif0t)";
            return;
        }
        for (int h = 0; h < myharmonics.size(); h++)
        {
            // Make sure the harmonic exists:
            if (myharmonics[h].size() != 0)
            {
                int freqindex = harmonic::getfrequency(h);
                // If the harmonic is a sine:
                if (harmonic::issine(h))
                    std::cout << " +vs" << freqindex << "*sin(" << freqindex*2 << "*pif0t)";
                else
                    std::cout << " +vc" << freqindex << "*cos(" << freqindex*2 << "*pif0t)";
            }
        }
        std::cout << std::endl;
    }
}

std::shared_ptr<coefmanager> rawfield::resetcoefmanager(void)
{
    synchronize();
    
    std::shared_ptr<coefmanager> outcm = mycoefmanager;
    
    mycoefmanager = std::shared_ptr<coefmanager>(new coefmanager(mytypename, myptracker->getdisjointregions()));
    
    for (int i = 0; i < interpolationorder.size(); i++)
        mycoefmanager->fitinterpolationorder(i, interpolationorder[i]);
    
    return outcm;
}

void rawfield::setcoefmanager(std::shared_ptr<coefmanager> cm)
{
    synchronize();
    
    mycoefmanager = cm;
}

void rawfield::print(void)
{
    // Do not sync this.
    
    if (myname.size() == 0)
        std::cout << "field";
    else
        std::cout << myname;
}

void rawfield::printvalues(bool databoundsonly)
{
    synchronize();
    
    if (mycoefmanager == NULL)
        std::cout << "Field has no values to print." << std::endl;
    else
        mycoefmanager->print(databoundsonly);
}

void rawfield::setname(std::string name)
{
    // Do not sync this.
    
    myname = name;
    
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setname(name);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setname(name);
    }
}

std::string rawfield::gettypename(bool familyonly)
{
    // Do not sync this.
    
    std::string out = mytypename;

    if (familyonly == false)
    {
        if (mysubfields.size() == 2)
            out = out+"xy";
        if (mysubfields.size() == 3)
            out = out+"xyz";
    }
    
    return out;
}

void rawfield::setorder(int physreg, int interpolorder, bool iscalledbyuser)
{
    if (iscalledbyuser && ispadaptive)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: .setorder(physreg, interpolorder) cannot be called on fields set to p-adaptivity" << std::endl;
        log.error();
    }
    
    // Interpolation order can only be set on highest dimension regions:
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    int regdim = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getelementdimension();
    if (iscalledbyuser && regdim >= 0 && problemdimension != regdim)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set the interpolation order on a " << regdim << "D region in a " << problemdimension << "D problem (must be " << problemdimension << "D)" << std::endl;
        log.error();
    }
    // Interpolation orders must be provided in decreasing order by the user to guarantee field continuity at the interfaces:
    if (iscalledbyuser && myordertracker.size() > 0 && myordertracker[myordertracker.size()-1].second < interpolorder)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: interpolation orders must be set descendingly to guarantee field continuity at the region interfaces" << std::endl;
        log.error();
    }

    synchronize();
    
    // Keep track of the calls to 'setorder':
    if (issynchronizing == false && mysubfields.size() == 0 && myharmonics.size() == 0)
        myordertracker.push_back(std::make_pair(physreg, interpolorder));
        
    if (mytypename == "x" || mytypename == "y" || mytypename == "z" || mytypename == "one0" || mytypename == "one1" || mytypename == "one2" || mytypename == "one3")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot choose the interpolation order for the x, y, z coordinate or for 'one' type fields" << std::endl;
        log.error();
    }
        
    // Set the interpolation order on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setorder(physreg, interpolorder, iscalledbyuser);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setorder(physreg, interpolorder, iscalledbyuser);
    }
        
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
    {
        // Consider ALL disjoint regions in the physical region with (-1):
        std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);
        
        hierarchicalformfunction myhff;
        int lowestfieldorder = myhff.getminorder(mytypename);
    
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            if (isitported[selecteddisjregs[i]])
            {
                interpolationorder[selecteddisjregs[i]] = lowestfieldorder;
                mycoefmanager->fitinterpolationorder(selecteddisjregs[i], lowestfieldorder);
            }
            else
            {
                interpolationorder[selecteddisjregs[i]] = interpolorder;
                mycoefmanager->fitinterpolationorder(selecteddisjregs[i], interpolorder);
            }
        }
    }
}

void rawfield::setorder(expression criterion, int loworder, int highorder, double critrange)
{
    synchronize();
    
    if (mytypename == "x" || mytypename == "y" || mytypename == "z" || mytypename == "one0" || mytypename == "one1" || mytypename == "one2" || mytypename == "one3")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot choose the interpolation order for the x, y, z coordinate or for 'one' type fields" << std::endl;
        log.error();
    }
    
    // Set the interpolation order on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setorder(criterion, loworder, highorder, critrange);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setorder(criterion, loworder, highorder, critrange);
    }

    if (mysubfields.size() == 0 && myharmonics.size() == 0)
    {
        ispadaptive = true;
        
        criterion = sl::abs(criterion);
        
        universe::getrawmesh()->add(shared_from_this(), criterion, loworder, highorder, critrange);
    }
}

void rawfield::setport(int physreg, std::shared_ptr<rawport> primal, std::shared_ptr<rawport> dual)
{
    synchronize();
 
    if (sl::isempty(physreg))
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set a port to empty physical region " << physreg << std::endl;
        log.error();
    }
    
    if (mytypename != "h1" && mytypename != "hcurl")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set ports to '" << mytypename << "' type fields" << std::endl;
        log.error();
    }
    if (mysubfields.size() > 0)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set ports to fields with multiple components (work with individual components instead)" << std::endl;
        log.error();
    }
    
    std::vector<int> fieldharms = getharmonics();
    if (fieldharms != primal->getharmonics() || fieldharms != dual->getharmonics())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set a port with a harmonic content that does not match the field" << std::endl;
        log.error();
    }
    
    hierarchicalformfunction myhff;
    int lowestfieldorder = myhff.getminorder(mytypename);
    
    std::vector<bool> isddmdisjreg = universe::getrawmesh()->getdtracker()->isddmdisjointregion();
    
    for (int h = 0; h < fieldharms.size(); h++)
    {
        int harm = fieldharms[h];
        std::shared_ptr<rawport> ph = primal->harmonic(harm);
        std::shared_ptr<rawport> dh = dual->harmonic(harm);
        std::shared_ptr<rawfield> fh = harmonic(harm);
    
        if (issynchronizing == false && (ph->isassociated() || dh->isassociated()))
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: at least one port to set is already associated to a field" << std::endl;
            log.error();
        }

        ph->associate(true, dh, physreg, fh);
        dh->associate(false, ph, physreg, fh);
        
        std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            int cdr = selecteddisjregs[i];
            
            if (isddmdisjreg[cdr])
            {
                logs log;
                log.msg() << "Error in 'rawfield' object: trying to set a port on a DDM interface (ports must be confined inside a DDM domain)" << std::endl;
                log.error();
            }

            // Set lowest order:
            fh->interpolationorder[cdr] = lowestfieldorder;
            fh->mycoefmanager->fitinterpolationorder(cdr, lowestfieldorder);
            
            fh->mydisjregconstraints[cdr] = NULL;
            fh->myconditionalconstraints[cdr] = {};
            fh->isitgauged[cdr] = false;
            fh->isitported[cdr] = true;
        }
        
        // Keep track of the calls to 'setport':
        if (issynchronizing == false)
            fh->myporttracker.push_back(std::make_tuple(physreg, ph, dh));
    }
}

void rawfield::setvalue(int physreg, int numfftharms, expression* meshdeform, expression input, int eio)
{
    synchronize();
    
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set the value for the x, y or z coordinate" << std::endl;
        log.error();
    }
    if (input.countcolumns() != 1 || input.countrows() != countcomponents())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: the rawfield value must be set with a " << countcomponents() << "x1 expression" << std::endl;
        log.error();
    }

    // Set the values on the subfields:
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->setvalue(physreg, numfftharms, meshdeform, input.at(i,0), eio);
        return;
    }
    // Set the values on the harmonics:
    if (myharmonics.size() > 0)
    {
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                myharmonics[h][0]->setvalue(physreg, numfftharms, meshdeform, sl::getharmonic(h, input, numfftharms), eio);
        }
        return;
    }
    
    physicalregion* curpr = getrawmesh()->getphysicalregions()->get(physreg);
    
    if (input.iszero())
        getpointer()->setzerovalue(physreg);
    else
        getpointer()->updateshapefunctions(input, meshdeform, {curpr->getdisjointregions(0), curpr->getdisjointregions(1), curpr->getdisjointregions(2), curpr->getdisjointregions(3)}, eio);
}

void rawfield::setvalue(int physreg)
{
    synchronize();
    
    switch (countcomponents())
    {
        case 1:
            setvalue(physreg, -1, NULL, 0);
            break;
        case 2:
            setvalue(physreg, -1, NULL, expression(2,1,{0,0}));
            break;
        case 3:
            setvalue(physreg, -1, NULL, expression(3,1,{0,0,0}));
            break;
    }
}

void rawfield::setvalue(elementselector& elemselect, std::vector<double>& gpcoordsin, expression* meshdeform, densemat values)
{
    synchronize();
    
    if (mytypename != "h1d0" && mytypename != "h1d1" && mytypename != "h1d2" && mytypename != "h1d3")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: expected a 'h1d' type field to set the value at given gauss points" << std::endl;
        log.error();
    }

    int numelems = elemselect.countinselection();
    int elementtypenumber = elemselect.getelementtypenumber();
    std::vector<int> elemnums = elemselect.getelementnumbers();
    
    // Needs to be monoharmonic without FFT (the field cannot store the time vals).
    if (values.countrows() != numelems || 3*values.countcolumns() != gpcoordsin.size())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: in 'setvalue' received a " << values.countrows() << "x" << values.countcolumns() << " densemat (expected " << numelems << "x" << gpcoordsin.size()/3 << ")" << std::endl;
        log.error();
    }
    
    gausspoints gp(elementtypenumber, gpcoordsin);

    int numgp = gp.count();
    std::vector<double> gpcoords = gp.getcoordinates();
    std::vector<double> gpweights = gp.getweights();
    
    // Manages reuse automatically:
    std::shared_ptr<opdetjac> dj(new opdetjac);
    densemat detjac = dj->interpolate(elemselect, gpcoords, meshdeform)[1][0];
    // The Jacobian determinant should be positive irrespective of the node numbering:
    detjac.abs();


    // Shape functions 'h1d' are not orientation dependent -> no orientation loop required.
    // Their order can however be non-unique in the selected elements.

    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();

    // Group disj. regs. with same interpolation order (all have same element type number).
    std::vector<int> interpolorders(alldisjregs.size());
    for (int i = 0; i < alldisjregs.size(); i++)
        interpolorders[i] = getinterpolationorder(alldisjregs[i]);
    disjointregionselector mydisjregselector(alldisjregs, {interpolorders});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);

        elemselect.selectdisjointregions(mydisjregs);
        
        int numinselection = elemselect.countinselection();
        if (numinselection == 0)
            continue;
            
        int fforder = getinterpolationorder(mydisjregs[0]);
        std::vector<int> elemindexes = elemselect.getelementindexes();
        
        densemat dj = detjac.extractrows(elemindexes);
        densemat vl = values.extractrows(elemindexes);
        

        ///// Projection formulation -> integral(dofv * tfv - vl * tfv) = 0

        vl.multiplyelementwise(dj);
        vl.transpose();
        
        // 1. Evaluate dofv and tfv * gpweights:
        
        hierarchicalformfunctioncontainer ffval = *(universe::gethff(mytypename, elementtypenumber, fforder, gpcoords));
        densemat testfunctionvalue = ffval.tomatrix(0, fforder, 0, 0); // total orientation is always 0 for h1d type
        densemat doffunctionvalue = testfunctionvalue.copy();
        
        testfunctionvalue.multiplycolumns(gpweights);
            
        // 2. Compute the tf * dof product [tfrow1*dofrow1 tfrow1*dofrow2 ... tfrow2*dofrow1 ...]:
        densemat testfuntimesdof = testfunctionvalue.multiplyallrows(doffunctionvalue);
        testfuntimesdof.transpose();

        // 3. Create matrix A terms:
        densemat Avals = dj.multiply(testfuntimesdof); // will be already row-col sorted
        
        // 4. Create right handside vector b terms:
        densemat bvals = testfunctionvalue.multiply(vl);

        // 5. Create A (block diagonal):
        int numffs = testfunctionvalue.countrows();
        int matsize = numinselection*numffs;
        
        if (numffs == 0)
            continue;
        
        indexmat csrrows(1, matsize+1, 0,numffs);
        indexmat csrcols(numinselection, numffs*numffs);
        int* captr = csrcols.getvalues();
        
        int index = 0;
        for (int e = 0; e < numinselection; e++)
        {
            for (int t = 0; t < numffs; t++)
            {
                for (int d = 0; d < numffs; d++)
                {
                    captr[index] = e*numffs+d;
                    index++;
                }
            }
        }
        
        Mat bdmat;
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, matsize, matsize, csrrows.getvalues(), csrcols.getvalues(), Avals.getvalues(), &bdmat);

        MatAssemblyBegin(bdmat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(bdmat, MAT_FINAL_ASSEMBLY);

        indexmat blocksizes(numinselection,1, numffs);
        densemat blockvals(numinselection*numffs, numffs);
        MatInvertVariableBlockDiagonal(bdmat, numinselection, blocksizes.getvalues(), blockvals.getvalues());
        
        MatDestroy(&bdmat);
        
        // 6. Solve block-diagonal system:
        densemat prod = blockvals.blockdiagonaltimesvector(blocksizes, bvals.gettranspose());
        double* prodptr = prod.getvalues();

        // 7. Set coefs to coefmanager:
        
        elements* els = universe::getrawmesh()->getelements();
        disjointregions* drs = universe::getrawmesh()->getdisjointregions();
        for (int e = 0; e < numinselection; e++)
        {
            int curel = elemnums[elemindexes[e]];
            int disjreg = els->getdisjointregion(elementtypenumber, curel);
            int rb = drs->getrangebegin(disjreg);
            int elemindexindr = curel-rb;

            for (int ff = 0; ff < numffs; ff++)
                mycoefmanager->setcoef(disjreg, ff, elemindexindr, prodptr[e*numffs+ff]);
        }
    }
    
    // Unselect the disjoint regions:
    elemselect.selectdisjointregions({});
}

void rawfield::setnodalvalues(indexmat nodenumbers, densemat values)
{
    synchronize();
    
    if (myharmonics.size() == 2 && myharmonics[1].size() == 1)
    {
        myharmonics[1][0]->setnodalvalues(nodenumbers, values);
        return;
    }
    
    elements* els = universe::getrawmesh()->getelements();
    disjointregions* drs = universe::getrawmesh()->getdisjointregions(); 

    int numnodes = nodenumbers.count();
    int* nptr = nodenumbers.getvalues();
    double* vptr = values.getvalues();

    // Only for 'h1' type fields:
    if (mytypename != "h1")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set nodal values for '" << mytypename << "' type fields (only 'h1' type)" << std::endl;
        log.error();
    }
    
    if (mysubfields.size() != 0 || myharmonics.size() != 0)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set nodal values for fields with subfields/harmonics (select a single component/harmonic)" << std::endl;
        log.error();
    }
    
    for (int i = 0; i < numnodes; i++)
    {
        int curnode = nptr[i];
        int ndr = els->getdisjointregion(0, curnode, false);

        if (ndr >= 0)
        {
            int rb = drs->getrangebegin(ndr);
            mycoefmanager->setcoef(ndr, 0, curnode-rb, vptr[i]);
        }
        else
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: cannot set value for node number " << curnode << " (is not an element corner node)" << std::endl;
            log.error();
        }
    }
}

densemat rawfield::getnodalvalues(indexmat nodenumbers)
{
    synchronize();
    
    if (myharmonics.size() == 2 && myharmonics[1].size() == 1)
        return myharmonics[1][0]->getnodalvalues(nodenumbers);
    
    elements* els = universe::getrawmesh()->getelements();
    disjointregions* drs = universe::getrawmesh()->getdisjointregions(); 
    
    int numnodes = nodenumbers.count();
    int* nptr = nodenumbers.getvalues();
    
    densemat output(nodenumbers.countrows(), nodenumbers.countcolumns());
    double* vptr = output.getvalues();

    // Only for 'h1' type fields:
    if (mytypename != "h1")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot get nodal values for '" << mytypename << "' type fields (only 'h1' type)" << std::endl;
        log.error();
    }
    
    if (mysubfields.size() != 0 || myharmonics.size() != 0)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot get nodal values for fields with subfields/harmonics (select a single component/harmonic)" << std::endl;
        log.error();
    }
    
    for (int i = 0; i < numnodes; i++)
    {
        int curnode = nptr[i];
        int ndr = els->getdisjointregion(0, curnode, false);

        if (ndr >= 0)
        {
            int rb = drs->getrangebegin(ndr);
            vptr[i] = mycoefmanager->getcoef(ndr, 0, curnode-rb);
        }
        else
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: cannot get value for node number " << curnode << " (is not an element corner node)" << std::endl;
            log.error();
        }
    }

    return output;
}

void rawfield::setzerovalue(int physreg)
{
    synchronize();

    // Get all disjoint regions in the physical region:
    std::vector<int> disjregs = universe::getrawmesh()->getphysicalregions()->get(physreg)->getdisjointregions(-1);

    for (int i = 0; i < disjregs.size(); i++)
    {
        int numelems = universe::getrawmesh()->getdisjointregions()->countelements(disjregs[i]);

        for (int ff = 0; ff < mycoefmanager->countformfunctions(disjregs[i]); ff++)
        {
            for (int elem = 0; elem < numelems; elem++)
                mycoefmanager->setcoef(disjregs[i], ff, elem, 0);
        }
    }
}

void rawfield::setdisjregconstraint(int physreg, int numfftharms, expression* meshdeform, expression input, int eio)
{
    synchronize();
    
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot constrain the x, y or z coordinate" << std::endl;
        log.error();
    }
    if (input.countcolumns() != 1 || input.countrows() != countcomponents())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: the rawfield must be constrained using a " << countcomponents() << "x1 expression" << std::endl;
        log.error();
    }
        
    // Set the constraints on the subfields:
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->setdisjregconstraint(physreg, numfftharms, meshdeform, input.at(i,0), eio);
        return;
    }
    // Set the constraints on the harmonics:
    if (myharmonics.size() > 0)
    {
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                myharmonics[h][0]->setdisjregconstraint(physreg, numfftharms, meshdeform, sl::getharmonic(h, input, numfftharms), eio);
        }
        return;
    }

    std::vector<expression> mdv = {};
    if (meshdeform != NULL)
        mdv = {*meshdeform};

    if (issynchronizing == false)
        mydisjregconstrainttracker.push_back(std::make_tuple(physreg, numfftharms, mdv, input, eio));

    // Find the highest tag:
    int maxtag = -1;
    for (int i = 0; i < mydisjregconstraints.size(); i++)
    {
        if (mydisjregconstraints[i] != NULL)
            maxtag = std::max(maxtag, std::get<5>(*mydisjregconstraints[i]));
    }
    std::shared_ptr<std::tuple<int,int,std::vector<expression>,expression,int,int>> ci(new std::tuple<int,int,std::vector<expression>,expression,int,int>{physreg, numfftharms, mdv, input, eio, maxtag+1});

    // Consider all disjoint regions in the physical region with (-1):
    std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);
    
    for (int i = 0; i < selecteddisjregs.size(); i++)
    {
        // Ports have priority over the disjreg constraints!
        if (isitported[selecteddisjregs[i]] == false)
        {
            mydisjregconstraints[selecteddisjregs[i]] = ci;
            myconditionalconstraints[selecteddisjregs[i]] = {};
            isitgauged[selecteddisjregs[i]] = false;
        }
    }
}

void rawfield::setdisjregconstraint(int physreg)
{
    synchronize();
    
    switch (countcomponents())
    {
        case 1:
            setdisjregconstraint(physreg, -1, NULL, 0);
            break;
        case 2:
            setdisjregconstraint(physreg, -1, NULL, expression(2,1,{0,0}));
            break;
        case 3:
            setdisjregconstraint(physreg, -1, NULL, expression(3,1,{0,0,0}));
            break;
    }
}

void rawfield::setconditionalconstraint(int physreg, expression condexpr, expression valexpr)
{
    synchronize();
    
    // Multiharmonic expressions are not allowed.
    std::vector<int> prdisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions();
    if (not(condexpr.isharmonicone(prdisjregs)) || not(valexpr.isharmonicone(prdisjregs)))
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set a conditional constraint with multiharmonic arguments" << std::endl;
        log.error();
    }
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot constrain the x, y or z coordinate" << std::endl;
        log.error();
    }
    if (valexpr.countcolumns() != 1 || valexpr.countrows() != countcomponents())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: the rawfield must be constrained using a " << countcomponents() << "x1 expression" << std::endl;
        log.error();
    }
    if (condexpr.countcolumns() != 1 || condexpr.countrows() != 1)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: expected a scalar condition for the conditional constraint" << std::endl;
        log.error();
    }
    
    // Set the conditional constraints on the subfields:
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->setconditionalconstraint(physreg, condexpr, valexpr.at(i,0));
        return;
    }
    if (myharmonics.size() == 2 && myharmonics[1].size() == 1)
    {
        myharmonics[1][0]->setconditionalconstraint(physreg, condexpr, valexpr);
        return;
    }
    if (myharmonics.size() != 0)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set conditional constraints for fields with harmonics (select a single harmonic)" << std::endl;
        log.error();
    }
    
    // Keep track of the calls to 'setconditionalconstraint':
    if (issynchronizing == false)
        myconditionalconstrainttracker.push_back(std::make_tuple(physreg, condexpr, valexpr));
        
    // Consider only the NODAL disjoint regions in the physical region with (0):
    std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(0);

    for (int i = 0; i < selecteddisjregs.size(); i++)
    {
        // Ports and disjreg constraints have priority over the conditional constraints!
        if (isitported[selecteddisjregs[i]] == false && mydisjregconstraints[selecteddisjregs[i]] == NULL)
            myconditionalconstraints[selecteddisjregs[i]] = {condexpr, valexpr};
    }
}

void rawfield::setgauge(int physreg)
{
    synchronize();
    
    // Keep track of the calls to 'setgauge':
    if (issynchronizing == false && mysubfields.size() == 0 && myharmonics.size() == 0)
        mygaugetracker.push_back(physreg);
        
    // Set the gauge on the subfields (if any):
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setgauge(physreg);
    // Set the gauge on the harmonics:
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setgauge(physreg);
    }
    
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
    {
        // Get ALL disjoint regions in the physical region (not only ders, to remove grad type form functions).
        std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);

        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            // Ports and disjreg constraints have priority over the gauge!
            if (isitported[selecteddisjregs[i]] == false && mydisjregconstraints[selecteddisjregs[i]] == NULL)
                isitgauged[selecteddisjregs[i]] = true;
        }
    }
}

void rawfield::setspanningtree(std::shared_ptr<rawspanningtree> spantree)
{
    synchronize();
    
    // Set the spanning tree on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setspanningtree(spantree);
    // Set the spanning tree on the harmonics:
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setspanningtree(spantree);
    }
    
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
        myspanningtree = spantree;
}

std::shared_ptr<rawspanningtree> rawfield::getspanningtree(void)
{
    synchronize();
    
    if (myspanningtree != NULL)
        return myspanningtree;
    else
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: spanning tree was not provided to rawfield" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::shared_ptr<rawfield> rawfield::getpointer(void)
{
    synchronize();
   
    return shared_from_this();
}

std::shared_ptr<rawmesh> rawfield::getrawmesh(void)
{
    synchronize();
    
    return myrawmesh;
}

std::shared_ptr<ptracker> rawfield::getptracker(void)
{
    synchronize();
    
    return myptracker;
}

std::shared_ptr<coefmanager> rawfield::getcoefmanager(void)
{
    synchronize();
    
    return mycoefmanager;
}

void rawfield::setdata(int physreg, vectorfieldselect myvec, std::string op)
{
    synchronize();
    
    // Extract the info from the vector with selected field:
    std::shared_ptr<rawfield> selectedrawfield = myvec.getrawfield();
    std::shared_ptr<rawvec> selectedvec = myvec.getrawvector();
        
    disjointregions* dr = selectedvec->getptracker()->getdisjointregions();
        
    // The raw fields must be of the same type:
    if (mytypename != selectedrawfield->mytypename)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: .setdata can only transfer data between fields of same type" << std::endl;
        log.error();
    }
    // The raw fields must have a same number of subfields:
    if (mysubfields.size() != selectedrawfield->mysubfields.size())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: .setdata can only transfer data from fields with same number of subfields" << std::endl;
        log.error();
    }
    
    // Get the data of every subfield:
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->setdata(physreg, vectorfieldselect(selectedvec, selectedrawfield->mysubfields[i][0]), op);
    }
    else
    {
        // The raw fields must include the same harmonic numbers:
        if (getharmonics() != selectedrawfield->getharmonics())
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: .setdata can only transfer data from fields with same harmonic numbers" << std::endl;
            log.error();
        }
        
        // Get the data of every harmonic:
        if (myharmonics.size() > 0)
        {
            for (int h = 0; h < myharmonics.size(); h++)
            {
                if (myharmonics[h].size() > 0)
                    myharmonics[h][0]->setdata(physreg, vectorfieldselect(selectedvec, selectedrawfield->harmonic(h)), op);
            }
        }
        else
        {
            // Extract the actual field from non-multiharmonic fields:
            if (selectedrawfield->amimultiharmonic == false)
                selectedrawfield = selectedrawfield->harmonic(1);
            
            // Get the data for a single field.
            
            // Get ALL disjoint regions in the physical region:
            std::vector<int> selecteddisjregs;
            if (physreg != -1)
                selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);
            else
            {
                std::shared_ptr<dofmanager> dofmngr = selectedvec->getdofmanager();
                dofmngr->selectfield(selectedrawfield);
                selecteddisjregs = dofmngr->getdisjointregionsofselectedfield();
            }

            for (int i = 0; i < selecteddisjregs.size(); i++)
            {
                int disjreg = selecteddisjregs[i];

                // In case the order of this raw field is higher than the order of 
                // the selected raw field we have to set to zero the higher orders.
                if (op == "set" && getinterpolationorder(disjreg) > selectedrawfield->getinterpolationorder(disjreg))
                {
                    // Decrease the order to forget the higher orders...
                    mycoefmanager->fitinterpolationorder(disjreg, selectedrawfield->getinterpolationorder(disjreg));
                    // ... then reset the previous order:
                    mycoefmanager->fitinterpolationorder(disjreg, getinterpolationorder(disjreg));
                }

                int elementtypenumber = dr->getelementtypenumber(disjreg);
                int elementdimension = dr->getelementdimension(disjreg);
                int numelem = dr->countelements(disjreg);

                std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, mytypename);
                // The interpolation order for this field and the selected fields might be different.
                int numformfunctionsperelement = std::min(myformfunction->count(getinterpolationorder(disjreg), elementdimension, 0), myformfunction->count(selectedrawfield->getinterpolationorder(disjreg), elementdimension, 0));

                for (int ff = 0; ff < numformfunctionsperelement; ff++)
                {
                    densemat values = selectedvec->getvalues(selectedrawfield, disjreg, ff);
                    double* vals = values.getvalues();
                    
                    // Transfer nothing if 'values' is empty:
                    if (vals != NULL)
                    {
                        if (op == "set")
                        {
                            for (int elem = 0; elem < numelem; elem++)
                                mycoefmanager->setcoef(disjreg, ff, elem, vals[elem]);
                        }
                        if (op == "add")
                        {
                            for (int elem = 0; elem < numelem; elem++)
                                mycoefmanager->setcoef(disjreg, ff, elem, mycoefmanager->getcoef(disjreg, ff, elem) + vals[elem]);
                        }
                    }
                }
            }
        }
    }
}

void rawfield::transferdata(int physreg, vectorfieldselect myvec, std::string op)
{
    synchronize();
    
    // Extract the info from the vector with selected field:
    std::shared_ptr<rawfield> selectedrawfield = myvec.getrawfield();
    std::shared_ptr<rawvec> selectedvec = myvec.getrawvector();
        
    // The raw fields must be of the same type:
    if (mytypename != selectedrawfield->mytypename)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: .transferdata can only transfer data between fields of same type" << std::endl;
        log.error();
    }
    // The raw fields must have a same number of subfields:
    if (mysubfields.size() != selectedrawfield->mysubfields.size())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: .transferdata can only transfer data from fields with same number of subfields" << std::endl;
        log.error();
    }
    
    // Transfer the data of every subfield:
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->transferdata(physreg, vectorfieldselect(selectedvec, selectedrawfield->mysubfields[i][0]), op);
        return;
    }
    
    // The raw fields must include the same harmonic numbers:
    if (getharmonics() != selectedrawfield->getharmonics())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: .transferdata can only transfer data from fields with same harmonic numbers" << std::endl;
        log.error();
    }
    
    // Transfer the data of every harmonic:
    if (myharmonics.size() > 0)
    {
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                myharmonics[h][0]->transferdata(physreg, vectorfieldselect(selectedvec, selectedrawfield->harmonic(h)), op);
        }
        return;
    }
    // Extract the actual field from non-multiharmonic fields:
    if (selectedrawfield->amimultiharmonic == false)
        selectedrawfield = selectedrawfield->harmonic(1);
    
    // Transfer the data for a single field.
    
    // Get ALL disjoint regions in the physical region:
    std::vector<int> selecteddisjregs;
    if (physreg != -1)
        selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);
    else
    {
        std::shared_ptr<dofmanager> dofmngr = selectedvec->getdofmanager();
        dofmngr->selectfield(selectedrawfield);
        selecteddisjregs = dofmngr->getdisjointregionsofselectedfield();
    }

    for (int i = 0; i < selecteddisjregs.size(); i++)
    {
        int disjreg = selecteddisjregs[i];

        int elementtypenumber = (universe::getrawmesh()->getdisjointregions())->getelementtypenumber(disjreg);
        int elementdimension = (universe::getrawmesh()->getdisjointregions())->getelementdimension(disjreg);
        int numelem = (universe::getrawmesh()->getdisjointregions())->countelements(disjreg);

        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, mytypename);
        // The interpolation order for this field and the selected fields might be different.
        int numformfunctionsinoriginfield = myformfunction->count(getinterpolationorder(disjreg), elementdimension, 0);
        int numformfunctionsperelement = myformfunction->count(selectedrawfield->getinterpolationorder(disjreg), elementdimension, 0);

        for (int ff = 0; ff < numformfunctionsperelement; ff++)
        {
            densemat values(1,numelem,0.0);
            double* vals = values.getvalues();
            
            if (ff < numformfunctionsinoriginfield)
            {
                for (int elem = 0; elem < numelem; elem++)
                    vals[elem] = mycoefmanager->getcoef(disjreg, ff, elem);
            }

            selectedvec->setvalues(selectedrawfield, disjreg, ff, values, op);
        }
    }
}

void rawfield::setcohomologysources(std::vector<int> cutphysregs, std::vector<double> cutvalues)
{
    synchronize();

    // There is no subfield for hcurl. No subfield loop needed.
    if (myharmonics.size() == 2 && myharmonics[1].size() > 0)
    {
        myharmonics[1][0]->setcohomologysources(cutphysregs, cutvalues);
        return;
    }
    
    if (myharmonics.size() > 0)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot set a cohomology source on a multiharmonic field (select harmonics one by one)" << std::endl;
        log.error();
    }
    
    elements* els = universe::getrawmesh()->getelements();
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    physicalregions* prs = universe::getrawmesh()->getphysicalregions();

    // Reset all cut sources:
    std::vector<std::vector<int>> ders(cutphysregs.size());
    for (int i = 0; i < cutphysregs.size(); i++)
    {
        ders[i] = prs->get(cutphysregs[i])->getdisjointregions(1);
        
        for (int d = 0; d < ders[i].size(); d++)
        {
            int ne = drs->countelements(ders[i][d]);
            int nff = mycoefmanager->countformfunctions(ders[i][d]);

            for (int f = 0; f < nff; f++)       
            {
                for (int j = 0; j < ne; j++)
                    mycoefmanager->setcoef(ders[i][d], f, j, 0.0);
            }   
        }
    }
    
    // Set all cut sources:
    std::vector<bool> flipit; std::vector<int> edgenums;
    for (int i = 0; i < cutphysregs.size(); i++)
    {
        if (ders[i].size() == 0)
            continue;
    
        gentools::inoutorient(cutphysregs[i], flipit);

        int index = 0;
        for (int d = 0; d < ders[i].size(); d++)
        {
            int rb = drs->getrangebegin(ders[i][d]);
            int ne = drs->countelements(ders[i][d]);
            
            for (int j = 0; j < ne; j++)
            {
                int curorient = 2 * els->gettotalorientation(1, rb+j) - 1;
                
                int signfix = curorient;
                if (flipit[index])
                    signfix *= -1;

                // Add the source with correct sign:
                double curval = mycoefmanager->getcoef(ders[i][d], 0, j);
                mycoefmanager->setcoef(ders[i][d], 0, j, curval + signfix * cutvalues[i]);
                
                index++;
            }
        }
    }
}

std::shared_ptr<rawfield> rawfield::comp(int component)
{   
    synchronize();
    
    // If there is a single component and the first one is requested:
    if (countcomponents() == 1 && component == 0)
        return shared_from_this();
    
    if (countformfunctioncomponents() > 1)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot get a component for vector fields with no subfields (e.g. hcurl)" << std::endl;
        log.error();
    }
    if (component > mysubfields.size())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot get component " << component << " from a " << mysubfields.size() << " components field" << std::endl;
        log.error();
    }
        
    return mysubfields[component][0];
}

std::shared_ptr<rawfield> rawfield::harmonic(int harmonicnumber)
{
    // Do not sync this.
    
    return harmonic(std::vector<int>{harmonicnumber});
}

std::shared_ptr<rawfield> rawfield::harmonic(const std::vector<int> harmonicnumbers)
{
    // Do not sync this.
    
    if (mysubfields.size() != 0)
    {
        std::shared_ptr<rawfield> harmsrawfield(new rawfield());
        *harmsrawfield = *this;
        // Replace every subfield by a new one with only the requested harmonics:
        for (int i = 0; i < mysubfields.size(); i++)
            harmsrawfield->mysubfields[i][0] = mysubfields[i][0]->harmonic(harmonicnumbers);
        return harmsrawfield;
    }
    if (myharmonics.size() != 0)
    {
        if (harmonicnumbers.size() == 1)
        {
            // If there is a single harmonic we can not put it into a new rawfield!
            if (harmonicnumbers[0] < myharmonics.size() && myharmonics[harmonicnumbers[0]].size() > 0)
                return myharmonics[harmonicnumbers[0]][0];
            else
            {
                logs log;
                log.msg() << "Error in 'rawfield' object: in .harmonic cannot get harmonic " << harmonicnumbers[0] << " (does not exist)" << std::endl; 
                log.error();
            }
        }
        
        // In case we want several harmonics a new raw field container is created.
        int maxharmnum = *max_element(harmonicnumbers.begin(), harmonicnumbers.end());

        std::shared_ptr<rawfield> harmsrawfield(new rawfield());
        *harmsrawfield = *this;

        // Set a brand new harmonic vector:
        harmsrawfield->myharmonics = std::vector<std::vector<std::shared_ptr<rawfield>>>(maxharmnum+1, std::vector<std::shared_ptr<rawfield>>(0));
        
        // Add only the requested harmonics:
        for (int i = 0; i < harmonicnumbers.size(); i++)
        {
            if (harmonicnumbers[i] < myharmonics.size() && myharmonics[harmonicnumbers[i]].size() > 0)
                harmsrawfield->myharmonics[harmonicnumbers[i]] = {myharmonics[harmonicnumbers[i]][0]};
            else
            {
                logs log;
                log.msg() << "Error in 'rawfield' object: in .harmonic cannot get harmonic " << harmonicnumbers[i] << " (does not exist)" << std::endl; 
                log.error();
            }
        }
        return harmsrawfield;
    }
        
    // In case there is no harmonic it is a constant (harmonic 1).
    if (harmonicnumbers.size() == 1 && harmonicnumbers[0] == 1)
        return shared_from_this();
    else
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: in .harmonic cannot get harmonic in constant field (does not exist)" << std::endl; 
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<std::shared_ptr<rawfield>> rawfield::getsons(void)
{
    // Do not sync this.

    if (mysubfields.size() > 0)
    {
        std::vector<std::shared_ptr<rawfield>> output = {};
        for (int i = 0; i < mysubfields.size(); i++)
        {
            std::vector<std::shared_ptr<rawfield>> curoutput = mysubfields[i][0]->getsons();
            output.insert(output.end(), curoutput.begin(), curoutput.end());
        }
        return output;
    }
    if (myharmonics.size() > 0)
    {
        std::vector<std::shared_ptr<rawfield>> output = {};
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                output.push_back(myharmonics[h][0]);
        }
        return output;
    }
    
    return {getpointer()};
}

bool rawfield::isdisjregconstrained(int disjreg)
{
    synchronize();
    
    return not(mydisjregconstraints[disjreg] == NULL);
}

std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>> rawfield::getdisjregconstraints(void)
{
    synchronize();
    
    return mydisjregconstraints;
}

bool rawfield::isconditionallyconstrained(int disjreg)
{
    synchronize();
    
    return (myconditionalconstraints[disjreg].size() > 0);
}

std::vector<std::vector<expression>> rawfield::getconditionalconstraints(void)
{
    synchronize();
    
    return myconditionalconstraints;
}

bool rawfield::isgauged(int disjreg) 
{ 
    synchronize();
    
    return isitgauged[disjreg];
}

int rawfield::getinterpolationorder(int disjreg) 
{ 
    synchronize();
    
    if (mytypename == "x" || mytypename == "y" || mytypename == "z" || mytypename == "one0" || mytypename == "one1" || mytypename == "one2" || mytypename == "one3")
        return 1;
        
    if (mysubfields.size() == 0)
    {
        int toreturn;
    
        if (myharmonics.size() == 0)
            toreturn = interpolationorder[disjreg];
        else
        {
            errornotsameinterpolationorder(disjreg);
            toreturn = myharmonics[getfirstharmonic()][0]->interpolationorder[disjreg];
        }
        if (toreturn == -1)
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: interpolation order is undefined on the region" << std::endl;
            log.msg() << "Define it with field.setorder(region, order)" << std::endl;
            log.error();
        }
        return toreturn;
    }
    else
    { 
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot get the interpolation order of a field with subfields" << std::endl;
        log.error();
    }
    
    throw std::runtime_error(""); // fix return warning
}

std::vector<int> rawfield::getinterpolationorders(void)
{
    synchronize();
    
    return interpolationorder;
}

int rawfield::getinterpolationorders(int elementtypenumber, std::vector<int>& elementnumbers, std::vector<int>& fieldorders)
{
    synchronize();
    
    int numelems = elementnumbers.size();
    fieldorders.resize(numelems);
    
    elements* els = universe::getrawmesh()->getelements();

    int maxorder = -1;
    for (int i = 0; i < numelems; i++)
    {
        int curelem = elementnumbers[i];
        int curdisjreg = els->getdisjointregion(elementtypenumber, curelem);
        
        int curorder = interpolationorder[curdisjreg];
        if (curorder > maxorder)
            maxorder = curorder;
        
        if (curorder != -1)
            fieldorders[i] = curorder;
        else
        {
            logs log;
            log.msg() << "Error in 'rawfield' object: interpolation order is undefined on the region" << std::endl;
            log.msg() << "Define it with field.setorder(region, order)" << std::endl;
            log.error();
        }
    }
    
    return maxorder;
}

void rawfield::getinterpolationorders(int fieldorder, double alpha, double absthres, std::vector<double>& weightsforeachorder, std::vector<int>& lowestorders)
{
    synchronize();
    
    int numorders = 1+fieldorder;
    int numelems = weightsforeachorder.size()/numorders;
    
    // Deduce the output:
    lowestorders = std::vector<int>(numelems);
    
    for (int i = 0; i < numelems; i++)
    {
        double totalweight = 0.0;
        for (int o = 0; o < numorders; o++)
            totalweight += weightsforeachorder[i*numorders+o];
            
        // Lowest order if total weight is too small:
        if (totalweight <= std::abs(absthres))
        {
            lowestorders[i] = 0; // min order is artificially brought down to 0 for h1 type as well (see weight calc)
            continue;
        }
            
        double weightthreshold = std::abs(alpha)*totalweight;
    
        double accumulatedweight = 0.0;
        for (int o = 0; o < numorders; o++)
        {
            lowestorders[i] = o;
        
            accumulatedweight += weightsforeachorder[i*numorders+o];

            if (accumulatedweight >= weightthreshold)
                break;
        }
    }
}

void rawfield::getweightsforeachorder(int elementtypenumber, int fieldorder, std::vector<int>& elementnumbers, std::vector<double>& weightsforeachorder)
{
    synchronize();
    
    int numelems = elementnumbers.size();
    int numorders = 1+fieldorder;
    
    weightsforeachorder = std::vector<double>(numelems*numorders, 0.0);

    std::vector<double> averagevals = {};
    if (mytypename == "h1" || mytypename == "h1d0" || mytypename == "h1d1" || mytypename == "h1d2" || mytypename == "h1d3")
    {
        getaverage(elementtypenumber, elementnumbers, 1, averagevals); // nodal shape functions are at order 1
        for (int i = 0; i < numelems; i++)
            weightsforeachorder[i*numorders+0] = std::abs(averagevals[i]);
    }

    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    elements* els = universe::getrawmesh()->getelements();
    
    element myelement(elementtypenumber);
    
    // Create a form function iterator to iterate through all form functions of the element type.
    hierarchicalformfunctioniterator myiterator(mytypename, elementtypenumber, fieldorder);

    for (int ff = 0; ff < myiterator.count(); ff++)
    {
        int associatedelementtype = myiterator.getassociatedelementtype();
        int formfunctionindex = myiterator.getformfunctionindexinnodeedgefacevolume();
        int formfunctionorder = myiterator.getformfunctionorder();
        
        int num = myiterator.getnodeedgefacevolumeindex();
        // For quad subelements in prisms and pyramids:
        if ((elementtypenumber == 6 || elementtypenumber == 7) && associatedelementtype == 3)
            num -= myelement.counttriangularfaces(); 

        for (int i = 0; i < numelems; i++)
        {
            int elem = elementnumbers[i];
            int currentsubelem = els->getsubelement(associatedelementtype, elementtypenumber, elem, num);
            int currentdisjointregion = els->getdisjointregion(associatedelementtype, currentsubelem);
            int rb = drs->getrangebegin(currentdisjointregion);
            
            double curcoef = mycoefmanager->getcoef(currentdisjointregion, formfunctionindex, currentsubelem-rb);
            if (formfunctionorder == 1 && averagevals.size() > 0)
                curcoef -= averagevals[i];
            
            weightsforeachorder[i*numorders+formfunctionorder] += std::abs(curcoef);
        }
        myiterator.next();
    }    
}

void rawfield::errornotsameinterpolationorder(int disjreg)
{
    synchronize();
    
    if (mysubfields.size() == 0)
    {
        if (myharmonics.size() > 0)
        {
            std::vector<int> harms = getharmonics();
            for (int h = 0; h < harms.size(); h++)
            {
                if (myharmonics[harms[h]][0]->getinterpolationorder(disjreg) == myharmonics[harms[0]][0]->getinterpolationorder(disjreg))
                    continue;
                else
                {
                    logs log;
                    log.msg() << "Error in 'rawfield' object: the interpolation order must be the same for all harmonics" << std::endl;
                    log.error();
                }
            }
        }
    }
    else
    { 
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot call 'errornotsameinterpolationorder' on a field with subfields" << std::endl;
        log.error();
    }
}

void rawfield::getaverage(int elementtypenumber, std::vector<int>& elementnumbers, int maxorder, std::vector<double>& averagevals)
{
    synchronize();

    int numelems = elementnumbers.size();
    averagevals = std::vector<double>(numelems, 0.0); // init all zero
    element myelement(elementtypenumber);
    
    disjointregions* drs = universe::getrawmesh()->getdisjointregions();
    elements* els = universe::getrawmesh()->getelements();
    
    // Create a form function iterator to iterate through the form functions of the element type.
    hierarchicalformfunctioniterator myiterator(mytypename, elementtypenumber, maxorder);
    double invnumffs = 1.0/myiterator.count();

    for (int ff = 0; ff < myiterator.count(); ff++)
    {
        int associatedelementtype = myiterator.getassociatedelementtype();
        int formfunctionindex = myiterator.getformfunctionindexinnodeedgefacevolume();
        
        int num = myiterator.getnodeedgefacevolumeindex();
        // For quad subelements in prisms and pyramids:
        if ((elementtypenumber == 6 || elementtypenumber == 7) && associatedelementtype == 3)
            num -= myelement.counttriangularfaces(); 

        for (int i = 0; i < numelems; i++)
        {
            int elem = elementnumbers[i];
            int currentsubelem = els->getsubelement(associatedelementtype, elementtypenumber, elem, num);
            int currentdisjointregion = els->getdisjointregion(associatedelementtype, currentsubelem);
            int rb = drs->getrangebegin(currentdisjointregion);
            
            double curcoef = mycoefmanager->getcoef(currentdisjointregion, formfunctionindex, currentsubelem-rb);
            
            averagevals[i] += curcoef*invnumffs;
        }
        myiterator.next();
    }
}


std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> rawfield::getallsons(void)
{
    synchronize();
    
    std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> output = {};
    
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
        {
            if (mysubfields[i].size() > 0)
            {
                std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> cur = mysubfields[i][0]->getallsons();
                for (int f = 0; f < cur.size(); f++)
                {
                    cur[f].first[0] = i;
                    output.push_back(cur[f]);
                }
            }
        }
        return output;
    }
    
    if (myharmonics.size() > 0)
    {
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
            {
                std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> cur = myharmonics[h][0]->getallsons();
                for (int f = 0; f < cur.size(); f++)
                {
                    cur[f].first[1] = h;
                    output.push_back(cur[f]);
                }
            }
        }
        return output;
    }
    
    std::vector<int> curfirst = {0,1};
    output = {std::make_pair(curfirst, getpointer())};
    
    return output;
}

void rawfield::writeraw(int physreg, std::string filename, bool isbinary, std::vector<double> extradata)
{
    synchronize();
    
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: cannot write field type '" << mytypename << "' to raw format" << std::endl;
        log.error();
    }
    
    int numsubfields = countsubfields();
    std::vector<int> harmoniclist = getharmonics();
    int numharms = harmoniclist.size();

    hierarchicalformfunction myhff;
    int fieldtypenum = myhff.gettypenumber(mytypename);

    // Get all disjoint regions in the physical region:
    std::vector<int> selecteddisjregs = ((universe::getrawmesh()->getphysicalregions())->get(physreg))->getdisjointregions(-1);
    int numdisjregs = selecteddisjregs.size();
    disjointregions* mydisjointregions = universe::getrawmesh()->getdisjointregions();

    // Get all sons:
    std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> allsons = getallsons();
    int numsons = allsons.size();
    
    
    // Preallocate the int data vector. It consists in:
    //
    // 1. File format number
    // 2. Field type number
    // 3. Number of subfields
    // 4. Number of harmonics in the field
    // 5. Harmonic numbers in the field (must be identical for every subfield)
    // 6. Number of disjoint regions
    // 7. All disjoint region numbers interlaced with the number of elements in each disjoint region
    // 8. For every son loop on every disjoint region and store {interpolorder, doublevecrangebegin}
    //
    std::vector<int> intdata(5 + numharms + 2*numdisjregs + numsons * 2*numdisjregs);

    intdata[0] = 1; // Format number
    intdata[1] = fieldtypenum;
    intdata[2] = numsubfields;
    intdata[3] = numharms;
    
    int indexinintvec = 4;
    for (int i = 0; i < numharms; i++)  
        intdata[indexinintvec+i] = harmoniclist[i];
    indexinintvec += numharms;
    
    intdata[indexinintvec] = numdisjregs;
    indexinintvec++;
    for (int i = 0; i < numdisjregs; i++)
    {
        int disjreg = selecteddisjregs[i];
        intdata[indexinintvec+0] = disjreg;
        intdata[indexinintvec+1] = mydisjointregions->countelements(disjreg);
        indexinintvec += 2;
    }
    
    // Temporary double data container.
    // One vector per son and per disjoint region (son listing order is first subfield then each harmonic of the first subfield, then next subfield...):
    std::vector<std::vector<std::vector<double>>> doubledata(numsons, std::vector<std::vector<double>>(numdisjregs));


    int doubledatasize = extradata.size();
    for (int i = 0; i < numsons; i++)
    {
        std::shared_ptr<rawfield> curson = allsons[i].second;
        std::shared_ptr<coefmanager> curcoefmanager = curson->mycoefmanager;
    
        // For every son provide info for each disjoint region as {order, doublevecrangebegin}:
        for (int d = 0; d < numdisjregs; d++)
        {
            int disjreg = selecteddisjregs[d];
        
            int interpolorder = curson->getinterpolationorder(disjreg);
         
            // Get the element type number, element dimension and number of elements in the current disjoint region:
            int elementtypenumber = mydisjointregions->getelementtypenumber(disjreg);
            int elementdimension = mydisjointregions->getelementdimension(disjreg);
            int numberofelements = mydisjointregions->countelements(disjreg);
            
            // Get the number of form functions associated to dimension elementdimension:
            std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, mytypename);
            int numberofformfunctions = myformfunction->count(interpolorder, elementdimension, 0);

            // Populate the int data vector:
            intdata[indexinintvec+0] = interpolorder; 
            intdata[indexinintvec+1] = doubledatasize;
            indexinintvec += 2;
            
            doubledatasize += numberofelements * numberofformfunctions;
            
            
            // Populate the 'doubledata' vector:
            doubledata[i][d].resize(numberofelements * numberofformfunctions);
                        
            int curindex = 0;
            for (int ff = 0; ff < numberofformfunctions; ff++)
            {
                for (int e = 0; e < numberofelements; e++)
                {
                    doubledata[i][d][curindex] = curcoefmanager->getcoef(disjreg, ff, e);
                    curindex++;
                }
            }
        }
    }
    
    
    // Flatten the double data vector:
    std::vector<double> flatdoubledata(doubledatasize);
    
    int curindex = 0;
    for (int i = 0; i < extradata.size(); i++)
    {
        flatdoubledata[curindex] = extradata[i];
        curindex++;
    }
    
    for (int i = 0; i < doubledata.size(); i++)
    {
        for (int j = 0; j < doubledata[i].size(); j++)
        {
            for (int k = 0; k < doubledata[i][j].size(); k++)
            {
                flatdoubledata[curindex] = doubledata[i][j][k];
                curindex++;
            }
        } 
    }
    
    
    // This will write the int data length and double data length at the file begin:
    iointerface::write(filename, intdata, flatdoubledata, isbinary);
}

std::vector<double> rawfield::loadraw(std::string filename, bool isbinary)
{
    synchronize();
    
    disjointregions* mydisjointregions = universe::getrawmesh()->getdisjointregions();
    
    
    ///// Load the data from disk:
    std::vector<int> intdata;
    std::vector<double> doubledata;

    iointerface::load(filename, intdata, doubledata, isbinary);
    
    
    ///// Extract the file format number:
    // int fileformatnum = intdata[0];
    
    
    ///// Extract the field type number and make sure the field type names match:
    int fieldtypenum = intdata[1];
    hierarchicalformfunction myhff;
    std::string fieldtypename = myhff.gettypename(fieldtypenum);

    if (mytypename != fieldtypename)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: trying to load a '" << fieldtypename << "' type field from file '" << filename << "' to a '" << mytypename << "' type" << std::endl;
        log.error();
    }
    
    
    ///// Extract the number of subfields and make sure they match:
    int numsubfields = intdata[2];
    if (numsubfields != countsubfields())
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: trying to load a field with " << numsubfields << " subfields from file '" << filename << "' to one with " << countsubfields() << " subfields" << std::endl;
        log.error();
    }
    
    
    ///// Extract the number of harmonics and all harmonics, make sure they match:
    int numharms = intdata[3];
    
    int indexinintvec = 4;
    std::vector<int> harmoniclist(numharms);
    for (int i = 0; i < numharms; i++)
        harmoniclist[i] = intdata[indexinintvec+i];
    indexinintvec += numharms;
    
    std::vector<int> harmsinthisfield = getharmonics();

    if (harmoniclist != harmsinthisfield)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: harmonic list in field loaded from file '" << filename << "' does not match current field." << std::endl << std::endl;

        log.msg() << "Harmonics in loaded field (" << harmoniclist.size() << "): ";
        for (int i = 0; i < harmoniclist.size(); i++)
            log.msg() << harmoniclist[i] << " ";
        log.msg() << std::endl;
        
        log.msg() << "Harmonics in current field (" << harmsinthisfield.size() << "): ";
        for (int i = 0; i < harmsinthisfield.size(); i++)
            log.msg() << harmsinthisfield[i] << " ";
        log.msg() << std::endl << std::endl;
        log.error();
    }
    
    
    ///// Extract the list of disjoint regions and the number of elements in them from intdata, make sure the meshes match:
    int numdisjregs = intdata[indexinintvec];
    indexinintvec++;
    std::vector<int> disjregs(numdisjregs); std::vector<int> numelems(numdisjregs);
    
    if (numdisjregs == 0)
        return doubledata;
        
    for (int i = 0; i < numdisjregs; i++)
    {
        disjregs[i] = intdata[indexinintvec+0];
        numelems[i] = intdata[indexinintvec+1];
        indexinintvec += 2;
    }
    int maxdisjreg = *max_element(disjregs.begin(), disjregs.end());
    
    bool issamemesh = true;
    if (maxdisjreg >= mydisjointregions->count())
        issamemesh = false;
    if (issamemesh == true)
    {
        for (int i = 0; i < numdisjregs; i++)
        {
            if (mydisjointregions->countelements(disjregs[i]) != numelems[i])
            {
                issamemesh = false;
                break;
            }
        }
    }
    if (issamemesh == false)
    {
        logs log;
        log.msg() << "Error in 'rawfield' object: the mesh used to write file '" << filename << "' is not the same as the current one" << std::endl;
        log.error();
    }
    
    
    ///// ALL CHECKS PERFORMED. LOAD THE DATA.
    
    // Load the extra data at the double data vector begin:
    std::vector<double> extradata(intdata[indexinintvec+1]);
    for (int i = 0; i < extradata.size(); i++)
        extradata[i] = doubledata[i];
    
    // Calculate the number of sons:
    int numsons = (intdata.size() - indexinintvec) / (2*numdisjregs);
    // Get all the sons in this rawfield:
    std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> allsons = getallsons();
    
    for (int i = 0; i < numsons; i++)
    {
        std::shared_ptr<rawfield> curson = allsons[i].second;
        std::shared_ptr<coefmanager> curcoefmanager = curson->mycoefmanager;
    
        for (int d = 0; d < numdisjregs; d++)
        {
            int disjreg = disjregs[d];
        
            int interpolorder = intdata[indexinintvec+0];
            int rangebeginindoubledata = intdata[indexinintvec+1];
            indexinintvec += 2;
            
            // Update the current field order:
            curcoefmanager->fitinterpolationorder(disjreg, interpolorder);
            curson->interpolationorder[disjreg] = interpolorder;
            
         
            // Get the element type number, element dimension and number of elements in the current disjoint region:
            int elementtypenumber = mydisjointregions->getelementtypenumber(disjreg);
            int elementdimension = mydisjointregions->getelementdimension(disjreg);
            int numberofelements = mydisjointregions->countelements(disjreg);
            
            // Get the number of form functions associated to dimension elementdimension:
            std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, mytypename);
            int numberofformfunctions = myformfunction->count(interpolorder, elementdimension, 0);

            
            // Extract the double data:
            int curindex = 0;
            for (int ff = 0; ff < numberofformfunctions; ff++)
            {
                for (int e = 0; e < numberofelements; e++)
                {
                    curcoefmanager->setcoef(disjreg, ff, e, doubledata[rangebeginindoubledata + curindex]);
                    curindex++;
                }
            }
        }
    }
    
    return extradata;
}


std::vector<densemat> rawfield::getjacterms(elementselector& elemselect, std::vector<double>& evaluationcoordinates)
{
    synchronize();
    
    int elementdimension = elemselect.getelementdimension();
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    int elementtypenumber = elemselect.getelementtypenumber();
    
    std::vector<densemat> output(elementdimension*problemdimension);

    elements* myelements = universe::getrawmesh()->getelements();

    std::vector<int> elementlist = elemselect.getelementnumbers();

    // Get all node coordinates:
    std::vector<double>* mynodecoordinates = universe::getrawmesh()->getnodes()->getcoordinates();

    element myelement(elementtypenumber, myelements->getcurvatureorder());        
    int numcurvednodes = myelement.countcurvednodes();

    int numrows = numcurvednodes;
    int numcols = elementlist.size();

    std::vector<densemat> coefmatrix(problemdimension);
    for (int d = 0; d < problemdimension; d++)
        coefmatrix[d] = densemat(numrows, numcols);
        
    std::vector<double*> coefs(problemdimension);
    for (int d = 0; d < problemdimension; d++)
        coefs[d] = coefmatrix[d].getvalues();

    for (int num = 0; num < numcurvednodes; num++)
    {
        for (int i = 0; i < elementlist.size(); i++)
        {
            int elem = elementlist[i];
            // Get the node to which the current form function is associated:
            int currentnode = myelements->getsubelement(0, elementtypenumber, elem, num);

            for (int d = 0; d < problemdimension; d++)
                coefs[d][num*numcols+i] = mynodecoordinates->at(3*currentnode+d);
        }
    }
    for (int d = 0; d < problemdimension; d++)
        coefmatrix[d].transpose();


    // Compute the form functions evaluated at the evaluation points:
    lagrangeformfunction mylagrange(elementtypenumber, myelements->getcurvatureorder(), evaluationcoordinates);
    std::vector<densemat> myformfunctionvalue(elementdimension);
    for (int ed = 0; ed < elementdimension; ed++)
        myformfunctionvalue[ed] = mylagrange.getderivative(1+ed);

    for (int ed = 0; ed < elementdimension; ed++)
    {
        for (int d = 0; d < problemdimension; d++)
            output[ed*problemdimension+d] = coefmatrix[d].multiply(myformfunctionvalue[ed]);
    }
    
    return output;
}
        

std::vector<std::vector<densemat>> rawfield::interpolate(int whichderivative, int formfunctioncomponent, elementselector& elemselect, std::vector<double>& evaluationcoordinates)
{
    synchronize();
    
    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();
    
    std::vector<std::vector<densemat>> out = {};

    // Group disj. regs. with same interpolation order (and same element type number).
    std::vector<int> interpolorders(alldisjregs.size());
    for (int i = 0; i < alldisjregs.size(); i++)
        interpolorders[i] = getinterpolationorder(alldisjregs[i]);
    disjointregionselector mydisjregselector(alldisjregs, {interpolorders});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
    
        elemselect.selectdisjointregions(mydisjregs);
        if (elemselect.countinselection() == 0)
            continue;
        
        int elementtypenumber = universe::getrawmesh()->getdisjointregions()->getelementtypenumber(mydisjregs[0]);
        std::vector<std::vector<densemat>> currentinterp = interpolate(whichderivative, formfunctioncomponent, elementtypenumber, elemselect.gettotalorientation(), getinterpolationorder(mydisjregs[0]), elemselect.getelementnumbers(), evaluationcoordinates);
        
        // Preallocate the harmonics not in 'out':
        if (out.size() < currentinterp.size())
            out.resize(currentinterp.size());
        for (int h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1 && out[h].size() == 0)
                out[h] = {densemat(elemselect.countincurrentorientation(), evaluationcoordinates.size()/3, 0)};
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
    
    return out;
}

densemat rawfield::getcoefficients(int elementtypenumber, int interpolorder, std::vector<int> elementlist)
{   
    synchronize();
    
    disjointregions* mydisjointregions = universe::getrawmesh()->getdisjointregions();  
    elements* myelements = universe::getrawmesh()->getelements();
    
    // In case the field is the x, y or z coordinate field:
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        int mycoordinate = 0;
        if (mytypename == "y")
            mycoordinate = 1;
        if (mytypename == "z")
            mycoordinate = 2;
        
        std::vector<double>* mynodecoordinates = universe::getrawmesh()->getnodes()->getcoordinates();
     
        element myelement(elementtypenumber, myelements->getcurvatureorder());        
        int numcurvednodes = myelement.countcurvednodes();
        
        int numrows = numcurvednodes;
        int numcols = elementlist.size();
        
        densemat coefmatrix(numrows, numcols);
        double* coefs = coefmatrix.getvalues();

        for (int num = 0; num < numcurvednodes; num++)
        {
            for (int i = 0; i < elementlist.size(); i++)
            {
                int elem = elementlist[i];
                // Get the node to which the current form function is associated:
                int currentnode = myelements->getsubelement(0, elementtypenumber, elem, num);
                
                coefs[num*numcols+i] = mynodecoordinates->at(3*currentnode+mycoordinate);
            }
        }
        return coefmatrix;
    }
    else
    {
        element myelement(elementtypenumber);
        
        // Create a form function iterator to iterate through all form functions.
        hierarchicalformfunctioniterator myiterator(mytypename, elementtypenumber, interpolorder);

        int numrows = myiterator.count();
        int numcols = elementlist.size();

        densemat coefmatrix(numrows, numcols);
        double* coefs = coefmatrix.getvalues();

        for (int ff = 0; ff < myiterator.count(); ff++)
        {
            int associatedelementtype = myiterator.getassociatedelementtype();
            int formfunctionindex = myiterator.getformfunctionindexinnodeedgefacevolume();
            int num = myiterator.getnodeedgefacevolumeindex();
            // For quad subelements in prisms and pyramids:
            if ((elementtypenumber == 6 || elementtypenumber == 7) && associatedelementtype == 3)
                num -= myelement.counttriangularfaces();

            for (int i = 0; i < elementlist.size(); i++)
            {
                int elem = elementlist[i];
                // Get the subelement to which the current form function is associated:
                int currentsubelem = myelements->getsubelement(associatedelementtype, elementtypenumber, elem, num);
                // Also get its disjoint region number:
                int currentdisjointregion = myelements->getdisjointregion(associatedelementtype, currentsubelem);
                // Use it to get the subelem index in the disjoint region:
                currentsubelem -= mydisjointregions->getrangebegin(currentdisjointregion);

                coefs[ff*numcols+i] = mycoefmanager->getcoef(currentdisjointregion, formfunctionindex, currentsubelem);
            }
            myiterator.next();
        }
        return coefmatrix;
    }
}

std::vector<std::vector<densemat>> rawfield::interpolate(int whichderivative, int formfunctioncomponent, int elementtypenumber, int totalorientation, int interpolorder, std::vector<int> elementnumbers, std::vector<double>& evaluationcoordinates)
{   
    synchronize();
    
    elements* myelements = universe::getrawmesh()->getelements();

    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {        
        // Get the coefficients. The interpolation order is not used here.
        densemat mycoefs = getcoefficients(elementtypenumber, -1, elementnumbers);
        // Compute the form functions evaluated at the evaluation points:
        lagrangeformfunction mylagrange(elementtypenumber, myelements->getcurvatureorder(), evaluationcoordinates);
        densemat myformfunctionvalue = mylagrange.getderivative(whichderivative);

        mycoefs.transpose();
        densemat coefstimesformfunctions = mycoefs.multiply(myformfunctionvalue);

        return {{},{coefstimesformfunctions}};
    }
    else
    {
        if (myharmonics.size() == 0)
        {
            // Get the coefficients:
            densemat mycoefs = getcoefficients(elementtypenumber, interpolorder, elementnumbers);
            // Compute the form functions evaluated at the evaluation points.
            // This reuses as much as possible what's already been computed:
            hierarchicalformfunctioncontainer* val = universe::gethff(mytypename, elementtypenumber, interpolorder, evaluationcoordinates);
            
            // We can get the total orientation from elementnumbers[0] since
            // we require that all elements have the same total orientation.
            densemat myformfunctionvalue = val->tomatrix(totalorientation, interpolorder, whichderivative, formfunctioncomponent);
            
            mycoefs.transpose();
            densemat coefstimesformfunctions = mycoefs.multiply(myformfunctionvalue);

            return {{},{coefstimesformfunctions}};
        }
        else
        {
            std::vector<std::vector<densemat>> coefstimesformfunctions(myharmonics.size());

            for (int harm = 1; harm < myharmonics.size(); harm++)
            {
                if (myharmonics[harm].size() == 0)
                    continue;
                else
                    coefstimesformfunctions[harm] = {(myharmonics[harm][0]->interpolate(whichderivative, formfunctioncomponent, elementtypenumber, totalorientation, interpolorder, elementnumbers, evaluationcoordinates))[1][0]};
            }
            return coefstimesformfunctions;
        }
    }
}


