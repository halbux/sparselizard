#include "dofinterpolate.h"


dofinterpolate::dofinterpolate(std::vector<double> refcoords, elementselector& elemselec, std::vector<std::shared_ptr<operation>> dofops, std::shared_ptr<dofmanager> dofmngr)
{
    oncontext* myoncontext = dofops[0]->getoncontext();

    mydoffield = dofops[0]->getfieldpointer();
    mydofops = dofops;
    mydofmanager = dofmngr;
    onphysreg = myoncontext->getphysicalregion();
    std::vector<int> ondisjregs = ((universe::getrawmesh()->getphysicalregions())->get(onphysreg))->getdisjointregions();
    myrefcoords = refcoords;
    mynumrefcoords = myrefcoords.size()/3;
    
    // Copy the element selector to keep it unchanged:
    elementselector elsel = elemselec;
    elsel.selectallelements();
    std::vector<int> originalindexes = elsel.getoriginalindexes();
    
    
    // Calculate the x, y and z coordinates of each reference coordinate to create the rcg object:
    myxyzcoords = std::vector<double>(3* elemselec.count() * mynumrefcoords);
    
    field x("x"), y("y"), z("z");
    expression xdef, ydef, zdef;
    xdef = x; ydef = y; zdef = z;
    
    if (myoncontext->isshifted())
    {
        expression shiftexpr = *(myoncontext->getshift());
    
        xdef = x+sl::compx(shiftexpr);
        if (shiftexpr.countrows() > 1)
            ydef = y+sl::compy(shiftexpr);
        if (shiftexpr.countrows() > 2)
            zdef = z+sl::compz(shiftexpr);
    }   

    bool isorientationdependent = (myoncontext->isshifted() && myoncontext->getshift()->isvalueorientationdependent(elsel.getdisjointregions()));
    std::vector<int> allelnums = elsel.getelementnumbers();
    elementselector cselsel(elsel.getdisjointregions(), allelnums, isorientationdependent);
    do
    {
        std::vector<int> csorigindexes = cselsel.getoriginalindexes();

        densemat xevaled = xdef.getoperationinarray(0,0)->interpolate(cselsel, myrefcoords, NULL)[1][0];
        densemat yevaled = ydef.getoperationinarray(0,0)->interpolate(cselsel, myrefcoords, NULL)[1][0];
        densemat zevaled = zdef.getoperationinarray(0,0)->interpolate(cselsel, myrefcoords, NULL)[1][0];

        double* xvals = xevaled.getvalues(); double* yvals = yevaled.getvalues(); double* zvals = zevaled.getvalues();

        int ind = 0;
        for (int e = 0; e < cselsel.countinselection(); e++)
        {
            for (int j = 0; j < mynumrefcoords; j++)
            {
                myxyzcoords[3*(originalindexes[csorigindexes[e]]*mynumrefcoords+j)+0] = xvals[ind];
                myxyzcoords[3*(originalindexes[csorigindexes[e]]*mynumrefcoords+j)+1] = yvals[ind];
                myxyzcoords[3*(originalindexes[csorigindexes[e]]*mynumrefcoords+j)+2] = zvals[ind];
                ind++;
            }
        }
    }
    while (cselsel.next());
        
    
    // Initialise to no coordinate found:
    isfound = std::vector<bool>(myxyzcoords.size()/3, false);

    rcg = referencecoordinategroup(myxyzcoords);


    // Calculate the maximum number of shape functions over all disjoint regions for preallocation.
    mymaxnumff = 0;
    for (int i = 0; i < ondisjregs.size(); i++)
    {
        int elementtypenumber = (universe::getrawmesh()->getdisjointregions())->getelementtypenumber(ondisjregs[i]); 
        int doforder = mydoffield->getinterpolationorder(ondisjregs[i]);
        std::shared_ptr<hierarchicalformfunction> dofformfunction = selector::select(elementtypenumber, mydoffield->gettypename());
        int curnumff = dofformfunction->count(doforder);
        if (mymaxnumff < curnumff)
            mymaxnumff = curnumff;
    }
    
    
    // Preallocate the matrix containers:
    int totalnumels = elemselec.count();
    myvals = std::vector<densemat>(mydofops.size());
    for (int i = 0; i < mydofops.size(); i++)
        myvals[i] = densemat(totalnumels, mymaxnumff*mynumrefcoords, 0);
    // Preallocate the adresses matrix for each harmonic:
    std::vector<int> dofharms = mydoffield->getharmonics();
    mydofnums = std::vector<std::vector<indexmat>>(*std::max_element(dofharms.begin(), dofharms.end()) + 1, std::vector<indexmat>(0));
   
    for (int h = 0; h < dofharms.size(); h++)
        mydofnums[dofharms[h]] = {indexmat(totalnumels, mymaxnumff*mynumrefcoords, -2)};
           

    eval();
}

void dofinterpolate::eval(void)
{
    elements* myelements = universe::getrawmesh()->getelements();
    disjointregions* mydisjointregions = universe::getrawmesh()->getdisjointregions();

    std::vector<int> ondisjregs = ((universe::getrawmesh()->getphysicalregions())->get(onphysreg))->getdisjointregions();

    std::vector<int> dofharms = mydoffield->getharmonics();
    
    std::vector<int> dofinterpolorders(ondisjregs.size());
    for (int i = 0; i < ondisjregs.size(); i++)
        dofinterpolorders[i] = mydoffield->getinterpolationorder(ondisjregs[i]);


    // Group disj. regs. with same element type and same dof interpolation order.
    disjointregionselector mydisjregselector(ondisjregs, {dofinterpolorders});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {    
        std::vector<int> disjregs = mydisjregselector.getgroup(i);
        
        // Calculate the reference coordinate positions:
        rcg.evalat(disjregs);
        
        int elementtypenumber = mydisjointregions->getelementtypenumber(disjregs[0]);        
        int doforder = mydoffield->getinterpolationorder(disjregs[0]);


        // Get all info on the shape functions:
        hierarchicalformfunctioniterator myiterator(mydoffield->gettypename(), elementtypenumber, doforder);
        int curnumff = myiterator.count();
        element myelement(elementtypenumber);
        std::vector<int> associatedelementtype(curnumff);
        std::vector<int> num(curnumff);
        std::vector<int> formfunctionindex(curnumff);
        for (int ff = 0; ff < curnumff; ff++)
        {
            associatedelementtype[ff] = myiterator.getassociatedelementtype();
            formfunctionindex[ff] = myiterator.getformfunctionindexinnodeedgefacevolume();
            num[ff] = myiterator.getnodeedgefacevolumeindex();
            // For quad subelements in prisms and pyramids:
            if ((elementtypenumber == 6 || elementtypenumber == 7) && associatedelementtype[ff] == 3)
                num[ff] -= myelement.counttriangularfaces();
                
            myiterator.next();
        }
        
        std::shared_ptr<hierarchicalformfunction> dofformfunction = selector::select(elementtypenumber, mydoffield->gettypename());
        hierarchicalformfunctioncontainer hfc = dofformfunction->evalat(doforder);

        bool isorientationdependent = dofformfunction->isorientationdependent(doforder);
        
        while (rcg.next())
        {
            std::vector<double> kietaphi = rcg.getreferencecoordinates();
            std::vector<int> coordindexes = rcg.getcoordinatenumber();
            std::vector<int> elemens = rcg.getelements();
            int numrefcoords = kietaphi.size()/3;
            
            // Evaluate the shape functions for all total orientations and all reference derivatives:
            hfc.evaluate(kietaphi);
            
            for (int h = 0; h < dofharms.size(); h++)
            {
                mydofmanager->selectfield(mydoffield->harmonic(dofharms[h]));
                int* dofnumsptr = mydofnums[dofharms[h]][0].getvalues();
                

                // Loop on all total orientations (if required):
                elementselector myselector(disjregs, elemens, isorientationdependent);
                do 
                {
                    int totalorient = myselector.gettotalorientation();
                    std::vector<int> elems = myselector.getelementnumbers();
                    std::vector<int> origindexes = myselector.getoriginalindexes();

                    // Rows are shape functions and columns are evaluation points:
                    std::vector<densemat> evaled(mydofops.size());
                    if (h == 0)
                    {
                        for (int df = 0; df < mydofops.size(); df++)
                            evaled[df] = hfc.tomatrix(totalorient, doforder, mydofops[df]->getkietaphiderivative(), mydofops[df]->getformfunctioncomponent());
                    }
                    
                    // Loop on all selected elements:
                    for (int e = 0; e < elems.size(); e++)
                    {
                        for (int ep = 0; ep < numrefcoords; ep++)
                        {
                            int callingindex = coordindexes[origindexes[e]*numrefcoords+ep];
                            
                            if (isfound[callingindex])
                                continue;
                            
                            int callingevalpt = callingindex%mynumrefcoords;
                            int callingelem = (callingindex-callingevalpt)/mynumrefcoords;
                            int rowstart = callingelem * mymaxnumff*mynumrefcoords;
   
                            // Same for all harmonics. Do it only once.
                            if (h == 0)
                            {
                                for (int df = 0; df < mydofops.size(); df++)
                                {
                                    for (int ff = 0; ff < curnumff; ff++)
                                        myvals[df].getvalues()[rowstart+mynumrefcoords*ff+callingevalpt] = evaled[df].getvalues()[ff*numrefcoords+ep];
                                }
                            }
                        
                            for (int ff = 0; ff < curnumff; ff++)
                            {
                                int currentsubelem = myelements->getsubelement(associatedelementtype[ff], elementtypenumber, elems[e], num[ff]);
                                int curdisjreg = myelements->getdisjointregion(associatedelementtype[ff], currentsubelem);
                                // In case of p-adaptivity some dofs can be missing and must be kept to the default negative adress.
                                if (mydofmanager->isdefined(curdisjreg, formfunctionindex[ff]))
                                {
                                    // Use it to get the subelem index in the disjoint region:
                                    currentsubelem -= mydisjointregions->getrangebegin(curdisjreg);
                                    
                                    if (mydofmanager->isported(curdisjreg))
                                        currentsubelem = 0;
                                    
                                    int rb = mydofmanager->getrangebegin(curdisjreg, formfunctionindex[ff]);
                                
                                    dofnumsptr[rowstart+mynumrefcoords*ff+callingevalpt] = rb + currentsubelem;
                                }
                            }
                        }
                    }         
                }
                while (myselector.next());
            }
            // Keep track of which coordinates were found:
            for (int c = 0; c < coordindexes.size(); c++)
                isfound[coordindexes[c]] = true;
        }
    }
    
    // Do something if some coordinates were not found:
    if (mydofops[0]->getoncontext()->iserrorifnotfound())
    {
        for (int i = 0; i < isfound.size(); i++)
        {    
            if (isfound[i] == false)
            {
                logs log;
                log.msg() << "Error in 'dofinterpolate' object: trying to interpolate at a point outside of region " << onphysreg << " or interpolation algorithm failed to converge" << std::endl;
                log.msg() << "Error was at (x,y,z) = (" << myxyzcoords[3*i+0] << ", " << myxyzcoords[3*i+1] << ", " << myxyzcoords[3*i+2] << ")" << std::endl;
                log.error();
            }
        }
    }
}


densemat dofinterpolate::getvalues(elementselector& elemselec, int dofopindex)
{
    std::vector<int> oi = elemselec.getoriginalindexes();
    return myvals[dofopindex].extractrows(oi);
}

indexmat dofinterpolate::getaddresses(elementselector& elemselec, int harmnum)
{
    std::vector<int> oi = elemselec.getoriginalindexes();
    return mydofnums[harmnum][0].extractrows(oi);
}


