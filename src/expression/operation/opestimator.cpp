#include "opestimator.h"


opestimator::opestimator(std::string estimatortype, std::shared_ptr<operation> arg)
{
    mytype = estimatortype;
    myarg = arg;
    
    if (mytype == "zienkiewiczzhu")
    {
        field estim("one");
        estim.noautomaticupdate();
        
        myvalue = std::shared_ptr<opfield>(new opfield(estim.harmonic(1).getpointer()));
    }
}

std::vector<std::vector<densemat>> opestimator::interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    if (meshdeform != NULL)
    {
        logs log;
        log.msg() << "Error in 'opestimator' object: estimator cannot be computed on a deformed mesh" << std::endl;
        log.error();
    }

    // Get the value from the universe if available and reuse is enabled:
    if (reuse && universe::isreuseallowed)
    {
        int precomputedindex = universe::getindexofprecomputedvalue(shared_from_this());
        if (precomputedindex >= 0) { return universe::getprecomputed(precomputedindex); }
    }
    
    bool wasreuseallowed = universe::isreuseallowed;
    universe::forbidreuse();
    
    // Update the estimator if allowed:
    if (universe::isestimatorupdateallowed(mystatenumber))
    {   
        if (mytype == "zienkiewiczzhu")
            estimatezienkiewiczzhu();
        mystatenumber = universe::estimatorcalcstate;
    }
    // If any update is forbidden reset value to 0:
    if (universe::numallowedtimes <= 0)
        myvalue->getfieldpointer()->resetcoefmanager();
    
    // Provide the requested output:
    std::vector<std::vector<densemat>> argmat = myvalue->interpolate(elemselect, evaluationcoordinates, NULL);
    
   if (wasreuseallowed)
        universe::allowreuse();
    
    if (reuse && universe::isreuseallowed)
        universe::setprecomputed(shared_from_this(), argmat);
    
    return argmat;
}

densemat opestimator::multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform)
{
    logs log;
    log.msg() << "Error in 'opestimator' object: cannot perform a multiharmonic interpolation on the estimator" << std::endl;
    log.error();
    
    throw std::runtime_error(""); // fix return warning
}

std::shared_ptr<operation> opestimator::simplify(std::vector<int> disjregs)
{
    myarg = myarg->simplify(disjregs);
    
    return shared_from_this();
}

std::shared_ptr<operation> opestimator::copy(void)
{
    std::shared_ptr<opestimator> op(new opestimator(mytype, myarg));
    // 'myvalue' should not be the same one as in this object thus do not do '*op = *this'.
    op->mystatenumber = 0; // to make sure 'myvalue' is recalculated when allowed (it is all zero here)
    return op;
}

void opestimator::print(void)
{
    std::cout << mytype << "(";
    myarg->print();
    std::cout << ")";
}

void opestimator::estimatezienkiewiczzhu(void)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    
    elements* myelements = universe::getrawmesh()->getelements();
    disjointregions* mydisjointregions = universe::getrawmesh()->getdisjointregions();
    
    // Compute 'myarg' at all mesh nodes: 
    std::vector<int> alldisjregsinmaxdim = mydisjointregions->getindim(problemdimension);

    int numnodes = universe::getrawmesh()->getnodes()->count();

    densemat nodalvaluesmin(numnodes, 1);
    densemat nodalvaluesmax(numnodes, 1);
    double* nvmin = nodalvaluesmin.getvalues();
    double* nvmax = nodalvaluesmax.getvalues();
    
    std::vector<int> numcontribs(numnodes, 0);

    // Send the disjoint regions with same element type numbers together:
    disjointregionselector mydisjregselector(alldisjregsinmaxdim, {});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);

        // Evaluate at the corner nodes:
        int elementtypenumber = mydisjointregions->getelementtypenumber(mydisjregs[0]);
    
        lagrangeformfunction mylagrange(elementtypenumber, 1, {});
        std::vector<double> evaluationpoints = mylagrange.getnodecoordinates();
        
        int nn = evaluationpoints.size()/3;

        // Loop on all total orientations (if required):
        bool isorientationdependent = myarg->isvalueorientationdependent(mydisjregs);
        elementselector myselector(mydisjregs, isorientationdependent);
        do
        {
            std::vector<int> elemnums = myselector.getelementnumbers();
        
            universe::allowreuse();
            densemat interpolated = myarg->interpolate(myselector, evaluationpoints, NULL)[1][0];
            universe::forbidreuse();

            double* interpvals = interpolated.getvalues();
            
            for (int e = 0; e < elemnums.size(); e++)
            {
                int curelem = elemnums[e];
                for (int n = 0; n < nn; n++)
                {
                    int curnode = myelements->getsubelement(0, elementtypenumber, curelem, n);
                    double curval = interpvals[e*nn+n];
                    
                    int numcont = numcontribs[curnode];
                
                    if (numcont == 0 || curval < nvmin[curnode])
                        nvmin[curnode] = curval;
                    if (numcont == 0 || curval > nvmax[curnode])
                        nvmax[curnode] = curval;
                    
                    numcontribs[curnode]++;
                }
            }
        }
        while (myselector.next());
    }
    
    // Populate 'myvalue':
    std::shared_ptr<rawfield> rf = myvalue->getfieldpointer();
    std::shared_ptr<coefmanager> cm = rf->getcoefmanager();
    
    for (int i = 0; i <= 7; i++)
    {
        int numelems = myelements->count(i);
        
        element elem(i);
        if (numelems == 0 || elem.getelementdimension() != problemdimension)
            continue;
    
        int nn = elem.countnodes();
    
        for (int e = 0; e < numelems; e++)
        {
            int curdisjreg = myelements->getdisjointregion(i, e);
            int rb = mydisjointregions->getrangebegin(curdisjreg);
        
            double curmaxerror = -1;
            for (int n = 0; n < nn; n++)
            {
                int curnode = myelements->getsubelement(0, i, e, n);
                double curerror = std::abs(nvmax[curnode]-nvmin[curnode]);
                if (curerror > curmaxerror)
                    curmaxerror = curerror;
            }
            
            cm->setcoef(curdisjreg, 0, e-rb, curmaxerror);
        }
    }
}
