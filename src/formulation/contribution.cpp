#include "contribution.h"


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

void contribution::generate(std::shared_ptr<rawvec> myvec, std::shared_ptr<rawmat> mymat, bool computeconstraints)
{   
    // Get the harmonics in the dof and tf fields. Get the max harmonic numbers as well.
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
        
    // The integration will be performed on the following disjoint regions:
    std::vector<int> selectedelemdisjregs = ((universe::mymesh->getphysicalregions())->get(integrationphysreg))->getdisjointregions();
  
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
    
    // Group disj. regs. with same element types and same tf and dof interpolation order.
    disjointregionselector mydisjregselector(selectedelemdisjregs, {tfinterpolorders, dofinterpolorders});
    for (int i = 0; i < mydisjregselector.countgroups(); i++)
    {
        std::vector<int> mydisjregs = mydisjregselector.getgroup(i);
        
        int elementtypenumber = (universe::mymesh->getdisjointregions())->getelementtypenumber(mydisjregs[0]);        
        
        // Get the interpolation order of the dof and tf fields:
        int tfinterpolationorder = tffield->getinterpolationorder(mydisjregs[0]);
        int dofinterpolationorder = 0;
        if (doffield != NULL)
            dofinterpolationorder = doffield->getinterpolationorder(mydisjregs[0]);
        
        // Compute the integration order, set it to zero if negative.
        // Adding an extra +2 generally gives a good integration in practice.
        int integrationorder = dofinterpolationorder + tfinterpolationorder + 2 + integrationorderdelta;
        if (integrationorder < 0)
        {
            std::cout << "Error in 'contribution' object: trying to integrate at negative order " << integrationorder << std::endl;
            abort();
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
        elementselector myselector(mydisjregs, isorientationdependent);
        do 
        {
            std::vector<int> elementnumbers = myselector.getelementnumbers();

            // stiffnesses[tf][dof][0] provides the stiffness matrix 
            // for test function harmonic 'tf' and dof harmonic 'dof'. 
            // If stiffnesses[tf][dof].size() is zero then it is empty.
            // stiffnesses[tf][1][0] must be used in case there is no dof.
            std::vector<std::vector<std::vector<densematrix>>> stiffnesses(maxtfharm + 1, std::vector<std::vector<densematrix>>(maxdofharm + 1, std::vector<densematrix>(0)));

            // Compute the Jacobian for the variable change to the reference element.
            std::shared_ptr<jacobian> myjacobian(new jacobian(myselector, evaluationpoints, meshdeformationptr));
            densematrix detjac = myjacobian->getdetjac();
            // The Jacobian determinant should be positive irrespective of the node numbering:
            detjac.abs();

            // Store it in the universe for reuse.
            universe::computedjacobian = myjacobian;
            universe::allowreuse();
            
            // Compute all terms in the contribution sum:
            densematrix tfformfunctionvalue, dofformfunctionvalue;
            for (int term = 0; term < mytfs.size(); term++)
            {
                ///// Compute the coefficients:
                // currentcoeff[i][0] holds the ith harmonic of the coefficient. 
                // It is empty if currentcoeff[i].size() is zero.
                std::vector<std::vector<densematrix>> currentcoeff;
                // Compute without or with FFT:
                if (numfftcoeffs <= 0)
                    currentcoeff = mycoeffs[term]->interpolate(myselector, evaluationpoints, meshdeformationptr);
                else
                {
                    densematrix timeevalinterpolated = mycoeffs[term]->multiharmonicinterpolate(numfftcoeffs, myselector, evaluationpoints, meshdeformationptr);
                    currentcoeff = myfft::fft(timeevalinterpolated, myselector.countinselection(), evaluationpoints.size()/3);
                }
                
                ///// Compute the dof*tf product (if any dof):
                densematrix doftimestestfun;
                tfformfunctionvalue = tfval.tomatrix(myselector.gettotalorientation(), tfinterpolationorder, mytfs[term]->getkietaphiderivative(), mytfs[term]->getformfunctioncomponent());
                        
                // Multiply by the weights:
                if (universe::skipgausspointweightproduct == false)
                    tfformfunctionvalue.multiplycolumns(weights);
                if (doffield != NULL)
                {
                    dofformfunctionvalue = dofval.tomatrix(myselector.gettotalorientation(), dofinterpolationorder, mydofs[term]->getkietaphiderivative(), mydofs[term]->getformfunctioncomponent());
                    doftimestestfun = tfformfunctionvalue.multiplyallrows(dofformfunctionvalue);
                }
                else
                    doftimestestfun = tfformfunctionvalue;

                ///// Since the interpolation orders are identical for all harmonics
                // we can premultiply all coefficients by the same dof*tf product.
                for (int h = 0; h < currentcoeff.size(); h++)
                {
                    if (currentcoeff[h].size() > 0)
                    {
                        currentcoeff[h][0].multiplyelementwise(detjac);
                        currentcoeff[h][0].transpose();
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
                                stiffnesses[currentharm][currentdofharm] = {currentcoeff[currentcoefharm][0].returnproduct(currentharmcoef)};
                            else
                                stiffnesses[currentharm][currentdofharm][0].addproduct(currentharmcoef, currentcoeff[currentcoefharm][0]);
                        }
                    }
                }
            }
            // Clear all reused data from the universe.
            universe::forbidreuse();
            
            // Get the adresses of all stiffnesses in the assembled matrix.
            for (int htf = 0; htf < tfharms.size(); htf++)
            {
                int currenttfharm = tfharms[htf];

                for (int hdof = 0; hdof < dofharms.size(); hdof++)
                {
                    int currentdofharm = dofharms[hdof];
                    
                    if (stiffnesses[currenttfharm][currentdofharm].size() == 0)
                        continue;
                        
                    ///// Get the adresses corresponding to every form function of 
                    // the test function/dof field in the elements of 'elementlist':
                    intdensematrix testfunadresses = mydofmanager->getadresses(tffield->harmonic(currenttfharm), tfinterpolationorder, elementtypenumber, elementnumbers, tfphysreg, computeconstraints);
                    intdensematrix dofadresses;
                    if (doffield != NULL)
                        dofadresses = mydofmanager->getadresses(doffield->harmonic(currentdofharm), dofinterpolationorder, elementtypenumber, elementnumbers, dofphysreg, false);

                    ///// Duplicate the tf and dof adresses to get an adress matrix of the size of the stiffness matrix.
                    if (doffield != NULL)
                    {
                        intdensematrix duplicateddofadresses = dofadresses.duplicateallrowstogether(tfformfunctionvalue.countrows());
                        intdensematrix duplicatedtestfunadresses = testfunadresses.duplicaterowsonebyone(dofformfunctionvalue.countrows());
                        
                        mymat->accumulate(duplicatedtestfunadresses, duplicateddofadresses, stiffnesses[currenttfharm][currentdofharm][0]);
                    }
                    else
                    {
                        // Bring back to the right hand side with a minus:
                        stiffnesses[currenttfharm][1][0].minus();
                        myvec->setvalues(testfunadresses, stiffnesses[currenttfharm][1][0], "add");
                        
                        // Keep track of how the rhs was assembled if requested:
                        if (universe::keeptrackofrhsassembly)
                            universe::rhsterms.push_back(std::make_pair(testfunadresses, stiffnesses[currenttfharm][1][0]));
                    }
                }
            }
        }
        while (myselector.next());        
    }
}
    
