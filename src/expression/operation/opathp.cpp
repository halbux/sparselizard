#include "opathp.h"
#include "myalgorithm.h"


opathp::opathp(std::shared_ptr<operation> arg, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt)
{
    myarg = arg;
    myrawmesh = rm;
    myptracker = pt;
}

std::vector<std::vector<densematrix>> opathp::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
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
                
    // All selected elements are of the same type:
    int numevalpts = evaluationcoordinates.size()/3;
    int elemtype = elemselect.getelementtypenumber();
    std::vector<int> elemnums = elemselect.getelementnumbers();
    
    int meshdim = universe::mymesh->getmeshdimension();
    
    
    ///// Get the corresponding reference coordinates on the highest dimension elements:
    
    std::vector<int> elemnumshere;
    std::vector<double> evalcoordshere;

    std::vector<int> maxdimdisjregs = universe::mymesh->getdisjointregions()->getindim(meshdim);
    universe::mymesh->getelements()->getrefcoordsondisjregs(elemtype, elemnums, evaluationcoordinates, maxdimdisjregs, elemnumshere, evalcoordshere);
    
    
    ///// Bring the evaluation points to the mesh here.
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


    ///// Bring the evaluation points to the ptracker here.
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
    
    
    ///// Evaluate the operation at all reference coordinate groups:

    densematrix argmat(elemnums.size(), numevalpts);
    double* argmatptr = argmat.getvalues();
    
    std::shared_ptr<rawmesh> bkp = universe::mymesh;
    universe::mymesh = myrawmesh->getattarget(myptracker);
    
    for (int i = 0; i < 8; i++)
    {
        element myelem(i);
        if (myelem.getelementdimension() != meshdim)
            continue;
    
        std::vector<int> curdisjregs = myptracker->getdisjointregions()->getintype(i);

        // Check if the operation is orientation dependent:
        bool isorientationdependent = myarg->isvalueorientationdependent(curdisjregs);
        
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
                densematrix interpoled = myarg->interpolate(myselector, kietaphi, meshdeform)[1][0];
                double* interpvals = interpoled.getvalues();
                      
                // Place the interpolated values at the right position in the output densematrix:
                std::vector<int> originds = myselector.getoriginalindexes();
                for (int j = 0; j < originds.size(); j++)
                {
                    int curorigelem = originds[j];
                    for (int k = 0; k < numrefcoords; k++)
                    {
                        int pos = coordindexes[curorigelem*numrefcoords+k];
                        argmatptr[pos] = interpvals[j*numrefcoords+k];
                    }
                }
            }
            while (myselector.next());  
        }
    }
    universe::mymesh = bkp;
    
    if (wasreuseallowed)
        universe::allowreuse();
    
    if (universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), {{}, {argmat}});
    
    return {{}, {argmat}};
}

void opathp::print(void)
{
    std::cout << "athp(";
    myarg->print();
    std::cout << ")";
}

