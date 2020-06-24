#include "opfieldnosync.h"
#include "myalgorithm.h"


opfieldnosync::opfieldnosync(std::shared_ptr<rawfield> fieldin)
{
    std::string tn = fieldin->gettypename();
    if (tn != "h1" && tn != "hcurl")
    {
        std::cout << "Error in 'opfieldnosync' object: cannot hp-adapt a '" << tn << "' type field" << std::endl;
        abort();
    }

    myfield = fieldin;
}

std::vector<std::vector<densematrix>> opfieldnosync::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    bool wasreuseallowed = universe::isreuseallowed;
    // Because of the 'gethff' call in interpolate:
    universe::forbidreuse();
    
    // Forbid synchronization:
    myfield->allowsynchronizing(false);
                
    // All selected elements are of the same type:
    int numevalpts = evaluationcoordinates.size()/3;
    int elemtype = elemselect.getelementtypenumber();
    int elemdim = elemselect.getelementdimension();
    std::vector<int> elemnums = elemselect.getelementnumbers();
    
    std::shared_ptr<rawmesh> myrawmesh = myfield->getrawmesh();
    std::shared_ptr<ptracker> myptracker = myfield->getptracker();
    
    
    ///// Bring the evaluation points to the not p-adapted universe::mymesh.
    //
    // UNIVERSE::MYMESH ---- P ----> UNIVERSE::MYMESH ---- h ----> myrawmesh ---- p ----> myptracker
    
    if (universe::mymesh != myrawmesh) 
    {
        std::vector<std::vector<int>> renumberingthere;
        universe::mymesh->getptracker()->getrenumbering(NULL, renumberingthere);

        for (int i = 0; i < elemnums.size(); i++)
            elemnums[i] = renumberingthere[elemtype][elemnums[i]];
    }

    
    ///// Bring the evaluation points to the mesh of this field.
    //
    // universe::mymesh ---- P ----> UNIVERSE::MYMESH ---- h ----> MYRAWMESH ---- p ----> myptracker
    
    std::vector<int> elemnumshere;
    std::vector<double> evalcoordshere;
    
    // In case there is no h-adaptivity:
    if (universe::mymesh == myrawmesh)
    {
        elemnumshere = std::vector<int>(2*elemnums.size()*numevalpts);
        for (int i = 0; i < elemnums.size(); i++)
        {
            // Give same format as when coming out of 'getattarget':
            for (int j = 0; j < numevalpts; j++)
            {
                elemnumshere[2*i*numevalpts+2*j+0] = elemtype;
                elemnumshere[2*i*numevalpts+2*j+1] = elemnums[i];
            }
        }
        evalcoordshere = myalgorithm::duplicate(evaluationcoordinates, elemnums.size());
    }
    else
    {
        std::vector<std::vector<int>> ad(8, std::vector<int>(0));
        std::vector<std::vector<double>> rc(8, std::vector<double>(0));
        
        int numelems = universe::mymesh->getelements()->count(elemtype);
        std::vector<bool> iselem(numelems, false);
        for (int i = 0; i < elemnums.size(); i++)
            iselem[elemnums[i]] = true;
        ad[elemtype] = std::vector<int>(numelems+1,0);
        for (int i = 1; i < numelems+1; i++)
        {
            if (iselem[i-1])
                ad[elemtype][i] = ad[elemtype][i-1] + evaluationcoordinates.size();
            else
                ad[elemtype][i] = ad[elemtype][i-1];
        }
        
        rc[elemtype] = myalgorithm::duplicate(evaluationcoordinates, elemnums.size());
        
        std::vector<std::vector<int>> tel;
        std::vector<std::vector<double>> trc;
        
        // Find at target:
        (universe::mymesh->gethtracker())->getattarget(ad, rc, myrawmesh->gethtracker().get(), tel, trc);
        
        // All element types have been placed in format type-number at position [elemtype]:
        elemnumshere = tel[elemtype];
        evalcoordshere = trc[elemtype];
    }

    
    ///// Bring the evaluation points to the ptracker of this field.
    //
    // universe::mymesh ---- P ----> universe::mymesh ---- h ----> MYRAWMESH ---- p ----> MYPTRACKER
    
    if (myrawmesh->getptracker() != myptracker)
    {
        std::vector<std::vector<int>> hererenumbering;
        (myrawmesh->getptracker())->getrenumbering(myptracker, hererenumbering);

        for (int i = 0; i < elemnumshere.size()/2; i++)
            elemnumshere[2*i+1] = hererenumbering[elemnumshere[2*i+0]][elemnumshere[2*i+1]];
    }


    ///// Place points in a 'referencecoordinategroup' object:
    
    referencecoordinategroup rcg(elemnumshere, evalcoordshere);
    
    
    ///// Evaluate the field at all reference coordinate groups:
    densematrix valmat(elemnums.size(), numevalpts);
    double* valsptr = valmat.getvalues();
    
    std::shared_ptr<rawmesh> bkp = universe::mymesh;
    universe::mymesh = myrawmesh->getattarget(myptracker);
    
    for (int i = 0; i < 8; i++)
    {
        element myelem(i);
        if (myelem.getelementdimension() != elemdim)
            continue;
        
        rcg.evalat(i);

        while (rcg.next()) //// SKIP ALL THIS IF ONLY P-ADAPTIVITY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            std::vector<double> kietaphi = rcg.getreferencecoordinates();
            std::vector<int> coordindexes = rcg.getcoordinatenumber();
            std::vector<int> elemens = rcg.getelements();
            int numrefcoords = kietaphi.size()/3;
            
            std::vector<int> curdisjregs = myptracker->getdisjointregions()->getintype(i);

            // Check if the field is orientation dependent:
            bool isorientationdependent = false;
            for (int j = 0; j < curdisjregs.size(); j++)
            {
                std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(j, myfield->gettypename());
                if ( myformfunction->isorientationdependent(myfield->getinterpolationorder(curdisjregs[j])) )
                    isorientationdependent = true;
            }

            // Loop on all total orientations (if required):
            elementselector myselector(curdisjregs, elemens, isorientationdependent);
            do 
            {
                densematrix interp = myfield->interpolate(0, formfunctioncomponent, myselector, kietaphi)[1][0];
                
                if (myfield->gettypename() == "hcurl")
                {
                    expression expr;
                    expression invjac = expr.invjac();
                    
                    densematrix fx = myfield->interpolate(0, 0, myselector, kietaphi)[1][0];
                    densematrix fy = myfield->interpolate(0, 1, myselector, kietaphi)[1][0];
                    densematrix fz = myfield->interpolate(0, 2, myselector, kietaphi)[1][0];
                    
                    densematrix invjacx = invjac.getoperationinarray(0,0)->interpolate(myselector, kietaphi, NULL)[1][0];
                    densematrix invjacy = invjac.getoperationinarray(0,1)->interpolate(myselector, kietaphi, NULL)[1][0];
                    densematrix invjacz = invjac.getoperationinarray(0,2)->interpolate(myselector, kietaphi, NULL)[1][0];
                    
                    fx.multiplyelementwise(invjacx);
                    fy.multiplyelementwise(invjacy);
                    fz.multiplyelementwise(invjacz);
                    
                    interp = densematrix(interp.countrows(), interp.countcolumns(), 0.0);
                    interp.add(fx);
                    interp.add(fy);
                    interp.add(fz);
                }
                
                double* interpvals = interp.getvalues();
                      
                // Place the interpolated values at the right position in the output densematrix:
                std::vector<int> originds = myselector.getoriginalindexes();
                for (int j = 0; j < originds.size(); j++)
                {
                    int curorigelem = originds[j];
                    for (int k = 0; k < numrefcoords; k++)
                    {
                        int pos = coordindexes[curorigelem*numrefcoords+k];
                        valsptr[pos] = interpvals[j*numrefcoords+k];
                    }
                }
            }
            while (myselector.next());  
        }
    }
    universe::mymesh = bkp;
    
    myfield->allowsynchronizing(true);
    
    std::vector<std::vector<densematrix>> out = {{}, {valmat}};
    
    if (wasreuseallowed)
        universe::allowreuse();
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), out);
    
    return out;
}

std::shared_ptr<operation> opfieldnosync::copy(void)
{
    std::shared_ptr<opfieldnosync> op(new opfieldnosync(myfield));
    *op = *this;
    return op;
}

void opfieldnosync::print(void)
{
    std::cout << "nosync(";
    myfield->print();
    std::cout << ")";
}

