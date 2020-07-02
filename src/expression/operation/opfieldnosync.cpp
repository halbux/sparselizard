#include "opfieldnosync.h"
#include "myalgorithm.h"


opfieldnosync::opfieldnosync(int formfunctioncomponent, std::shared_ptr<rawfield> fieldin)
{
    mycomp = formfunctioncomponent;
    myfield = fieldin;
}

void opfieldnosync::setcomponents(std::vector<std::shared_ptr<opfieldnosync>> allcomps)
{
    std::vector<std::weak_ptr<opfieldnosync>> weakptrs(allcomps.size());
    for (int c = 0; c < allcomps.size(); c++)
        weakptrs[c] = allcomps[c];

    mycomponents = weakptrs;
}

std::vector<std::vector<densematrix>> opfieldnosync::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{   
    // Get the value from the universe if available:
    if (universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }

    bool wasreuseallowed = universe::isreuseallowed;
    // Because of the 'gethff' call in interpolate:
    universe::forbidreuse();
    
    // Forbid synchronization:
    myfield->allowsynchronizing(false);
                
    std::string fieldtypename = myfield->gettypename();
                
    // All selected elements are of the same type:
    int numevalpts = evaluationcoordinates.size()/3;
    int elemtype = elemselect.getelementtypenumber();
    std::vector<int> elemnums = elemselect.getelementnumbers();
    
    int meshdim = universe::mymesh->getmeshdimension();
    
    std::shared_ptr<rawmesh> myrawmesh = myfield->getrawmesh();
    std::shared_ptr<ptracker> myptracker = myfield->getptracker();
    
    
    ///// Get the corresponding reference coordinates on the highest dimension elements:
    
    std::vector<int> elemnumshere;
    std::vector<double> evalcoordshere;

    std::vector<int> maxdimdisjregs = universe::mymesh->getdisjointregions()->getindim(meshdim);
    universe::mymesh->getelements()->getrefcoordsondisjregs(elemtype, elemnums, evaluationcoordinates, maxdimdisjregs, elemnumshere, evalcoordshere);
    
    
    ///// Bring the evaluation points to the mesh of this field.
    //
    // UNIVERSE::MYMESH ---- h ----> MYRAWMESH ---- p ----> myptracker
    
    if (universe::mymesh != myrawmesh)
    {
        std::vector<std::vector<int>> ads;
        std::vector<std::vector<double>> rcs;
        std::vector<int> indexinrcsoforigin;
        myalgorithm::toaddressdata(elemnumshere, evalcoordshere, universe::mymesh->getelements()->count(), ads, rcs, indexinrcsoforigin);
        
        // Find at target:
        std::vector<std::vector<int>> tel;
        std::vector<std::vector<double>> trc;
        (universe::mymesh->gethtracker())->getattarget(ads, rcs, myrawmesh->gethtracker().get(), tel, trc);
        
        // Recombine:
        for (int i = 0; i < elemnumshere.size()/2; i++)
        {
            int typ = elemnumshere[2*i+0];
            int ind = indexinrcsoforigin[i];
            
            elemnumshere[2*i+0] = tel[typ][2*ind+0];
            elemnumshere[2*i+1] = tel[typ][2*ind+1];

            evalcoordshere[3*i+0] = trc[typ][3*ind+0];
            evalcoordshere[3*i+1] = trc[typ][3*ind+1];
            evalcoordshere[3*i+2] = trc[typ][3*ind+2];
        }   
    }


    ///// Bring the evaluation points to the ptracker of this field.
    //
    // universe::mymesh ---- h ----> MYRAWMESH ---- p ----> MYPTRACKER
    
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
    
    int numcomps = mycomponents.size();
    std::vector<densematrix> valmats(numcomps);
    std::vector<double*> valsptrs(numcomps);
    for (int i = 0; i < numcomps; i++)
    {
        valmats[i] = densematrix(elemnums.size(), numevalpts);
        valsptrs[i] = valmats[i].getvalues();
    }
    
    std::shared_ptr<rawmesh> bkp = universe::mymesh;
    universe::mymesh = myrawmesh->getattarget(myptracker);
    
    for (int i = 0; i < 8; i++)
    {
        element myelem(i);
        if (myelem.getelementdimension() != meshdim)
            continue;
    
        std::vector<int> curdisjregs = myptracker->getdisjointregions()->getintype(i);

        // Check if the field is orientation dependent:
        bool isorientationdependent = false;
        for (int j = 0; j < curdisjregs.size(); j++)
        {
            std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(i, fieldtypename);
            if ( myformfunction->isorientationdependent(myfield->getinterpolationorder(curdisjregs[j])) )
                isorientationdependent = true;
        }
        
        rcg.evalat(i);

        while (rcg.next())
        {
            std::vector<double> kietaphi = rcg.getreferencecoordinates();
            std::vector<int> coordindexes = rcg.getcoordinatenumber();
            std::vector<int> elemens = rcg.getelements();
            int numrefcoords = kietaphi.size()/3;

            // Loop on all total orientations (if required):
            elementselector myselector(curdisjregs, elemens, isorientationdependent);
            do 
            {
                std::vector<densematrix> interpoled(numcomps);
                
                if (fieldtypename == "h1")
                    interpoled[0] = myfield->interpolate(0, 0, myselector, kietaphi)[1][0];
                
                if (fieldtypename == "hcurl")
                {
                    densematrix fx = myfield->interpolate(0, 0, myselector, kietaphi)[1][0];
                    densematrix fy = myfield->interpolate(0, 1, myselector, kietaphi)[1][0];
                    densematrix fz = myfield->interpolate(0, 2, myselector, kietaphi)[1][0];
                 
                    expression expr;
                    expression invjac = expr.invjac();
                    // To compute only once the jac:
                    universe::allowreuse();
                    for (int c = 0; c < numcomps; c++)
                    {
                        interpoled[c] = densematrix(myselector.countinselection(), kietaphi.size()/3, 0.0);
                        interpoled[c].addproduct( invjac.getoperationinarray(c, 0)->interpolate(myselector, kietaphi, NULL)[1][0], fx );
                        interpoled[c].addproduct( invjac.getoperationinarray(c, 1)->interpolate(myselector, kietaphi, NULL)[1][0], fy );
                        interpoled[c].addproduct( invjac.getoperationinarray(c, 2)->interpolate(myselector, kietaphi, NULL)[1][0], fz );
                    }
                    universe::forbidreuse();
                }
                
                std::vector<double*> interpvals(numcomps);
                for (int c = 0; c < numcomps; c++)
                    interpvals[c] = interpoled[c].getvalues();
                      
                // Place the interpolated values at the right position in the output densematrix:
                std::vector<int> originds = myselector.getoriginalindexes();
                for (int j = 0; j < originds.size(); j++)
                {
                    int curorigelem = originds[j];
                    for (int k = 0; k < numrefcoords; k++)
                    {
                        int pos = coordindexes[curorigelem*numrefcoords+k];
                        for (int c = 0; c < numcomps; c++)
                            valsptrs[c][pos] = interpvals[c][j*numrefcoords+k];
                    }
                }
            }
            while (myselector.next());  
        }
    }
    universe::mymesh = bkp;
    
    myfield->allowsynchronizing(true);
    
    if (wasreuseallowed)
        universe::allowreuse();
    
    if (universe::isreuseallowed)
    {
        for (int c = 0; c < numcomps; c++)
            universe::setprecomputed(mycomponents[c].lock(), {{}, {valmats[c]}});
    }
    
    return {{}, {valmats[mycomp]}};
}

void opfieldnosync::print(void)
{
    std::cout << "nosync(";
    myfield->print();
    std::cout << ")";
}

