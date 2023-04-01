#include "contribution.h"
#include "dofinterpolate.h"


contribution::contribution(std::shared_ptr<dofmanager> dofmngr) { mydofmanager = dofmngr; }

void contribution::setdofs(std::vector<std::shared_ptr<operation>> dofs) { mydofs = dofs; }
void contribution::settfs(std::vector<std::shared_ptr<operation>> tfs) { mytfs = tfs; }
void contribution::setcoeffs(std::vector<std::shared_ptr<operation>> coeffs) { mycoeffs = coeffs; }
void contribution::setdoffield(std::shared_ptr<rawfield> input) { doffield = input; }
void contribution::settffield(std::shared_ptr<rawfield> input) { tffield = input; }
void contribution::setmeshdeformation(expression meshdeform) { mymeshdeformation = {meshdeform}; }
void contribution::setintegrationphysicalregion(int physreg) { integrationphysreg = physreg; }
void contribution::setdofphysicalregion(int physreg) { dofphysreg = physreg; }
void contribution::settfphysicalregion(int physreg) { tfphysreg = physreg; }
void contribution::setintegrationorderdelta(int integrorderdelta) { integrationorderdelta = integrorderdelta; }
void contribution::setnumfftcoeffs(int numcoeffs) { numfftcoeffs = numcoeffs; }
void contribution::setbarycenterevalflag(void) { isbarycentereval = true; }

void contribution::generate(std::shared_ptr<rawvec> myvec, std::shared_ptr<rawmat> mymat)
{   
    bool isdofinterpolate = (doffield != NULL && mydofs[0]->ison());

    // Get the harmonics in the dof and tf fields. Get the max harmonic numbers as well:
    std::vector<int> tfharms = tffield->getharmonics();
    int maxtfharm = *std::max_element(tfharms.begin(), tfharms.end());
    std::vector<int> dofharms = {1};
    int maxdofharm = 1;
    if (doffield != NULL)
    {
        dofharms = doffield->getharmonics();
        maxdofharm = *std::max_element(dofharms.begin(), dofharms.end());   
    }
    // Get a pointer for the mesh deformation expression:
    expression* meshdeformationptr = NULL;
    if (mymeshdeformation.size() == 1)
        meshdeformationptr = &(mymeshdeformation[0]);
        
    // Check which node is in the tf physical region:
    std::vector<bool> isnodeintfphysreg;
    if (tfphysreg != integrationphysreg)
        universe::getrawmesh()->getelements()->istypeindisjointregions(0, universe::getrawmesh()->getphysicalregions()->get(tfphysreg)->getdefinition(), isnodeintfphysreg, false);
        
    // The integration will be performed on the following disjoint regions:
    std::vector<int> selectedelemdisjregs = universe::getrawmesh()->getphysicalregions()->get(integrationphysreg)->getdisjointregions();
  
    // Prepare to send the disjoint regions with same element type 
    // numbers and same dof and tf interpolation order together:
    std::vector<int> tfinterpolorders(selectedelemdisjregs.size());
    std::vector<int> dofinterpolorders(selectedelemdisjregs.size(),0);
    // All harmonics must have the same interpolation order!
    for (int i = 0; i < selectedelemdisjregs.size(); i++)
        tfinterpolorders[i] = tffield->getinterpolationorder(selectedelemdisjregs[i]);
    if (doffield != NULL)
    {
        for (int i = 0; i < selectedelemdisjregs.size(); i++)
            dofinterpolorders[i] = doffield->getinterpolationorder(selectedelemdisjregs[i]);
    }
    
    // Group disj. regs. with same element types and same tf and dof interpolation order:
    disjointregionselector mydisjregselector(selectedelemdisjregs, {tfinterpolorders, dofinterpolorders});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
        
        int elementtypenumber = (universe::getrawmesh()->getdisjointregions())->getelementtypenumber(mydisjregs[0]);        
        
        // Get the interpolation order of the dof and tf fields:
        int tfinterpolationorder = tffield->getinterpolationorder(mydisjregs[0]);
        int dofinterpolationorder = tfinterpolationorder;
        if (doffield != NULL)
            dofinterpolationorder = doffield->getinterpolationorder(mydisjregs[0]);
        
        // Compute the integration order.
        // Adding an extra +2 generally gives a good integration.
        int integrationorder = dofinterpolationorder + tfinterpolationorder + 2 + integrationorderdelta;
        if (isbarycentereval)
            integrationorder = 0;
        if (integrationorder < 0)
        {
            logs log;
            log.msg() << "Error in 'contribution' object: trying to integrate at negative order " << integrationorder << std::endl;
            log.error();
        }
            
        // Get the Gauss points and their weight:
        gausspoints mygausspoints(elementtypenumber, integrationorder);
        std::vector<double> evaluationpoints = mygausspoints.getcoordinates();
        std::vector<double> weights = mygausspoints.getweights();        
            
        // Compute the dof and tf form functions evaluated at the evaluation points:
        std::shared_ptr<hierarchicalformfunction> tfformfunction = selector::select(elementtypenumber, tffield->gettypename());
        
        // Copy needed here otherwise the pointed values will be overwritten by the dof:
        hierarchicalformfunctioncontainer tfval = *(universe::gethff(tffield->gettypename(), elementtypenumber, tfinterpolationorder, evaluationpoints));
        
        std::shared_ptr<hierarchicalformfunction> dofformfunction;
        hierarchicalformfunctioncontainer dofval;
        if (doffield != NULL)
        {
            dofformfunction = selector::select(elementtypenumber, doffield->gettypename());
            if (not(isdofinterpolate))
                dofval = *(universe::gethff(doffield->gettypename(), elementtypenumber, dofinterpolationorder, evaluationpoints));
        }
        
        // Simplify all coeffs for faster computation later on.
        // Also check if orientation matters.
        bool isorientationdependent = tfformfunction->isorientationdependent(tfinterpolationorder);
        if (doffield != NULL)
            isorientationdependent = isorientationdependent || dofformfunction->isorientationdependent(dofinterpolationorder);
        for (int term = 0; term < mytfs.size(); term++)
        {
            mycoeffs[term] = mycoeffs[term]->simplify(mydisjregs);
            isorientationdependent = (isorientationdependent || mycoeffs[term]->isvalueorientationdependent(mydisjregs) || (meshdeformationptr != NULL && meshdeformationptr->isvalueorientationdependent(mydisjregs)));
        }
        
        // Loop on all total orientations (if required):
        elementselector myselector;
        if (tfphysreg == integrationphysreg)
            myselector = elementselector(mydisjregs, isorientationdependent);
        else
        {
            std::vector<int> activeelems = gentools::getactiveelements(mydisjregs, isnodeintfphysreg);
            if (activeelems.size() > 0)
                myselector = elementselector(mydisjregs, activeelems, isorientationdependent);
        }
        if (myselector.count() == 0)
            continue;
        
        dofinterpolate mydofinterp;
        if (isdofinterpolate)
            mydofinterp = dofinterpolate(evaluationpoints, myselector, mydofs, mydofmanager);
        do 
        {
            std::vector<int> elementnumbers = myselector.getelementnumbers();

            // stiffnesses[tf][dof][0] provides the stiffness matrix 
            // for test function harmonic 'tf' and dof harmonic 'dof'. 
            // If stiffnesses[tf][dof].size() is zero then it is empty.
            // stiffnesses[tf][1][0] must be used in case there is no dof.
            std::vector<std::vector<std::vector<densemat>>> stiffnesses(maxtfharm + 1, std::vector<std::vector<densemat>>(maxdofharm + 1, std::vector<densemat>(0)));

            // Compute the Jacobian for the variable change to the reference element:
            std::shared_ptr<jacobian> myjacobian(new jacobian(myselector, evaluationpoints, meshdeformationptr));
            densemat detjac = myjacobian->getdetjac();
            // The Jacobian determinant should be positive irrespective of the node numbering:
            detjac.abs();

            // Store it in the universe for reuse:
            universe::computedjacobian = myjacobian;
            universe::allowreuse();
            
            // Compute all terms in the contribution sum:
            densemat tfformfunctionvalue, dofformfunctionvalue;
            for (int term = 0; term < mytfs.size(); term++)
            {
                ///// Compute the coefficients:
                // currentcoeff[i][0] holds the ith harmonic of the coefficient. 
                // It is empty if currentcoeff[i].size() is zero.
                std::vector<std::vector<densemat>> currentcoeff;
                // Compute without or with FFT:
                if (numfftcoeffs <= 0)
                    currentcoeff = mycoeffs[term]->interpolate(myselector, evaluationpoints, meshdeformationptr);
                else
                {
                    densemat timeevalinterpolated = mycoeffs[term]->multiharmonicinterpolate(numfftcoeffs, myselector, evaluationpoints, meshdeformationptr);
                    currentcoeff = fourier::fft(timeevalinterpolated, myselector.countinselection(), evaluationpoints.size()/3);
                }
                
                ///// Compute the dof*tf product (if any dof):
                densemat doftimestestfun;
                tfformfunctionvalue = tfval.tomatrix(myselector.gettotalorientation(), tfinterpolationorder, mytfs[term]->getkietaphiderivative(), mytfs[term]->getformfunctioncomponent());
                
                // Multiply by the weights:
                if (not(isbarycentereval))
                    tfformfunctionvalue.multiplycolumns(weights);
                if (doffield != NULL)
                {
                    if (isdofinterpolate)
                    {
                        dofformfunctionvalue = mydofinterp.getvalues(myselector, term);
                        doftimestestfun = dofformfunctionvalue.dofinterpoltimestf(tfformfunctionvalue);
                    }
                    else
                    {
                        dofformfunctionvalue = dofval.tomatrix(myselector.gettotalorientation(), dofinterpolationorder, mydofs[term]->getkietaphiderivative(), mydofs[term]->getformfunctioncomponent());
                        doftimestestfun = tfformfunctionvalue.multiplyallrows(dofformfunctionvalue);
                    }
                }
                else
                    doftimestestfun = tfformfunctionvalue;

                ///// Since the interpolation orders are identical for all harmonics
                // we can premultiply all coefficients by the same dof*tf product.
                for (int h = 0; h < currentcoeff.size(); h++)
                {
                    if (currentcoeff[h].size() > 0)
                    {
                        if (not(isbarycentereval))
                            currentcoeff[h][0].multiplyelementwise(detjac);
                        currentcoeff[h][0].transpose();
                        
                        if (isdofinterpolate)
                        {
                            densemat dttf = doftimestestfun.copy();
                            dttf.multiplycolumns(currentcoeff[h][0]);
                            currentcoeff[h][0] = dttf;  
                        }
                        else
                            currentcoeff[h][0] = doftimestestfun.multiply(currentcoeff[h][0]);
                    }
                }
                
                ///// Check if there is a time derivative on a multiharmonic dof:
                int multiharmonicdoftimederivativeorder = 0;
                if (doffield != NULL && doffield->ismultiharmonic())
                    multiharmonicdoftimederivativeorder = mydofs[term]->gettimederivative();

                ///// Add the term to the corresponding stiffness block:
                for (int currentcoefharm = 0; currentcoefharm < currentcoeff.size(); currentcoefharm++)
                {
                    if (currentcoeff[currentcoefharm].size() == 0)
                        continue;
                    
                    for (int dofharmindex = 0; dofharmindex < dofharms.size(); dofharmindex++)
                    {
                        int currentdofharm = dofharms[dofharmindex];
                        // Perform the product of the coefficient and the dof harmonic:
                        std::vector<std::pair<int, double>> harmsofproduct = harmonic::getproduct(currentcoefharm, currentdofharm, multiharmonicdoftimederivativeorder);
                        // Loop on all product harmonics:
                        for (int p = 0; p < harmsofproduct.size(); p++)
                        {
                            int currentharm = harmsofproduct[p].first;
                            // currentharmcoef can be + or - 0.5 or 1 (+ the time derivation factor).
                            double currentharmcoef = harmsofproduct[p].second;
                            
                            // Skip if the product harmonic is not a tf harmonic:
                            if (tffield->isharmonicincluded(currentharm) == false)
                                continue;
                            
                            // Add the term to the stiffnesses:
                            if (stiffnesses[currentharm][currentdofharm].size() == 0)
                                stiffnesses[currentharm][currentdofharm] = {currentcoeff[currentcoefharm][0].getproduct(currentharmcoef)};
                            else
                                stiffnesses[currentharm][currentdofharm][0].addproduct(currentharmcoef, currentcoeff[currentcoefharm][0]);
                        }
                    }
                }
            }
            // Clear all reused data from the universe:
            universe::forbidreuse();
            
            // Get the addresses of all stiffnesses in the assembled matrix:
            for (int htf = 0; htf < tfharms.size(); htf++)
            {
                int currenttfharm = tfharms[htf];

                for (int hdof = 0; hdof < dofharms.size(); hdof++)
                {
                    int currentdofharm = dofharms[hdof];
                    
                    if (stiffnesses[currenttfharm][currentdofharm].size() == 0)
                        continue;
                        
                    ///// Get the addresses corresponding to every form function of 
                    // the test function/dof field in the elements of 'elementlist':
                    indexmat testfunaddresses = mydofmanager->getaddresses(tffield->harmonic(currenttfharm), tfinterpolationorder, elementtypenumber, elementnumbers, tfphysreg);
                    indexmat dofaddresses;
                    if (doffield != NULL)
                    {
                        if (isdofinterpolate)
                            dofaddresses = mydofinterp.getaddresses(myselector, currentdofharm);
                        else
                            dofaddresses = mydofmanager->getaddresses(doffield->harmonic(currentdofharm), dofinterpolationorder, elementtypenumber, elementnumbers, dofphysreg);
                    }
                    
                    ///// Duplicate the tf and dof addresses as needed by the rawmat object:
                    if (doffield != NULL)
                    {
                        indexmat duplicateddofaddresses = dofaddresses;
                        if (isdofinterpolate)
                            duplicateddofaddresses = dofaddresses.duplicateallcolstogether(tfformfunctionvalue.countrows());
                            
                        indexmat duplicatedtestfunaddresses = testfunaddresses;
                        if (isdofinterpolate)
                        {
                            duplicatedtestfunaddresses = testfunaddresses.gettranspose();
                            duplicatedtestfunaddresses = duplicatedtestfunaddresses.duplicatecolsonebyone(dofaddresses.countcolumns());
                        }
                        
                        mymat->accumulate(duplicatedtestfunaddresses, duplicateddofaddresses, stiffnesses[currenttfharm][currentdofharm][0]);
                    }
                    else
                    {
                        // Bring back to the right hand side with a minus:
                        stiffnesses[currenttfharm][1][0].minus();
                        myvec->setvalues(testfunaddresses, stiffnesses[currenttfharm][1][0], "add");
                        
                        // Keep track of how the rhs was assembled if requested:
                        if (universe::keeptrackofrhsassembly)
                            universe::rhsterms.push_back(std::make_pair(testfunaddresses, stiffnesses[currenttfharm][1][0]));
                    }
                }
            }
        }
        while (myselector.next());        
    }
}
    
