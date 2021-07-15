#include "formulation.h"


formulation::formulation(void) { mydofmanager = std::shared_ptr<dofmanager>(new dofmanager); }

void formulation::operator+=(expression expr)
{
    if (isstructurelocked)
    {
        std::cout << "Error in 'formulation' object: cannot add port relations after a generation step" << std::endl;
        abort();
    }
    if (mycontributions[0].size() != 0 || mycontributions[1].size() != 0 || mycontributions[2].size() != 0 || mycontributions[3].size() != 0)
    {
        std::cout << "Error in 'formulation' object: port relations must be added before integral terms" << std::endl;
        abort();
    }

    // Add the port relation:
    std::shared_ptr<portrelation> pr(new portrelation(expr));
    myportrelations.push_back(pr);

    // Add all ports to the dofmanager:
    std::vector<std::shared_ptr<rawport>> rps = pr->getrawports();
    for (int p = 0; p < rps.size(); p++)
        mydofmanager->addtostructure(rps[p]);
}

void formulation::operator+=(std::vector<integration> integrationobject)
{
    for (int i = 0; i < integrationobject.size(); i++)
        *this += integrationobject[i];
}

void formulation::operator+=(integration integrationobject)
{
    if (isstructurelocked)
    {
        std::cout << "Error in 'formulation' object: cannot add contributions after a generation step" << std::endl;
        abort();
    }

    int integrationphysreg = integrationobject.getphysicalregion();
    int elementdimension = universe::mymesh->getphysicalregions()->get(integrationphysreg)->getelementdimension();
    // Return on empty integration region:
    if (elementdimension < 0)
        return;
    
    int integrationorderdelta = integrationobject.getintegrationorderdelta();
    expression myexpression = integrationobject.getexpression();

    myexpression.expand();
    
    // The element dimension is required to decompose space derivatives 
    // into the ki, eta and phi derivatives in the reference element.
    std::vector< std::vector<std::vector<std::shared_ptr<operation>>> > coeffdoftf = myexpression.extractdoftf(elementdimension);

    std::vector<std::vector<std::shared_ptr<operation>>> coeffs = coeffdoftf[0]; 
    std::vector<std::vector<std::shared_ptr<operation>>> dofs = coeffdoftf[1];
    std::vector<std::vector<std::shared_ptr<operation>>> tfs = coeffdoftf[2];

    // Loop on all slices:
    for (int slice = 0; slice < tfs.size(); slice++)
    {
        // In a given slice all dof and tf fields are the same, have the 
        // same applied time derivatives and are selected on a same 
        // physical region for all entries.
        std::shared_ptr<rawfield> doffield = dofs[slice][0]->getfieldpointer();
        std::shared_ptr<rawfield> tffield = tfs[slice][0]->getfieldpointer();
        
        // Get the time derivative of the dof field (if any) to add the
        // contribution to the K, C or M matrix. For a multiharmonic dof
        // the contribution is added to K.
        int contribindex = 0;
        if (doffield != NULL)
        {
            contribindex = 1;
            if (doffield->ismultiharmonic() == false)
                contribindex = dofs[slice][0]->gettimederivative() + 1;
        }
        
        // Get the physical regions on which the dof and tf are defined.
        int dofphysreg = -1;
        if (doffield != NULL)
            dofphysreg = dofs[slice][0]->getphysicalregion();
        int tfphysreg = tfs[slice][0]->getphysicalregion();
        // In case the physreg is -1 it is defined on the integration region:
        if (dofphysreg == -1)
            dofphysreg = integrationphysreg;
        if (tfphysreg == -1)
        tfphysreg = integrationphysreg;

        // Add the dofs to the dofmanager:
        if (doffield != NULL && dofs[slice][0]->ison() == false)
        {
            std::vector<int> dofharms = doffield->getharmonics();
            for (int h = 0; h < dofharms.size(); h++)
                mydofmanager->addtostructure(doffield->harmonic(dofharms[h]), dofphysreg);
        }
        std::vector<int> tfharms = tffield->getharmonics();
        for (int h = 0; h < tfharms.size(); h++)
            mydofmanager->addtostructure(tffield->harmonic(tfharms[h]), tfphysreg);

        // Create the contribution:
        contribution mycontribution(mydofmanager);

        mycontribution.setdoffield(doffield);
        mycontribution.settffield(tffield);

        mycontribution.setdofs(dofs[slice]);
        mycontribution.settfs(tfs[slice]);
        mycontribution.setcoeffs(coeffs[slice]);
        
        mycontribution.setintegrationorderdelta(integrationorderdelta);
        mycontribution.setnumfftcoeffs(integrationobject.getnumberofcoefharms());
        
        if (integrationobject.isbarycentereval)
            mycontribution.setbarycenterevalflag();
        
        mycontribution.setintegrationphysicalregion(integrationphysreg);
        if (doffield != NULL)
            mycontribution.setdofphysicalregion(dofphysreg);
        mycontribution.settfphysicalregion(tfphysreg);
        
        if (integrationobject.ismeshdeformdefined())
            mycontribution.setmeshdeformation(integrationobject.getmeshdeform());

        // Add the contribution to the contribution container:    
        int blocknumber = integrationobject.getblocknumber();
        if (mycontributions[contribindex].size() < blocknumber+1)
            mycontributions[contribindex].resize(blocknumber+1);
        mycontributions[contribindex][blocknumber].push_back(mycontribution);
    }
}

int formulation::countdofs(void)
{
    return mydofmanager->countdofs(); 
}

long long int formulation::allcountdofs(void)
{
    return mydofmanager->allcountdofs();
}

bool formulation::isstiffnessmatrixdefined(void)
{
    return (mycontributions[1].size() != 0);
}

bool formulation::isdampingmatrixdefined(void)
{
    return (mycontributions[2].size() != 0);
}

bool formulation::ismassmatrixdefined(void)
{
    return (mycontributions[3].size() != 0);
}

void formulation::generate(int m, int contributionnumber)
{
    isstructurelocked = true;
    
    if (contributionnumber < 0)
    {
        std::cout << "Error in 'formulation' object: cannot generate a negative contribution number" << std::endl;
        abort();
    }
    
    // Make sure the contribution number exists:
    if (contributionnumber >= mycontributions[m].size() || mycontributions[m][contributionnumber].size() == 0)
        return;
 
    universe::allowestimatorupdate(true);
        
    if (m == 0 && myvec == NULL)
        myvec = std::shared_ptr<rawvec>(new rawvec(mydofmanager));
    if (m > 0 && mymat[m-1] == NULL)
        mymat[m-1] = std::shared_ptr<rawmat>(new rawmat(mydofmanager));

    std::vector<contribution> contributionstogenerate = mycontributions[m][contributionnumber];
    for (int i = 0; i < contributionstogenerate.size(); i++)
    {
        if (m == 0)
            contributionstogenerate[i].generate(myvec, NULL);
        else
            contributionstogenerate[i].generate(NULL, mymat[m-1]);
    }
    
    universe::allowestimatorupdate(false);
    
}

void formulation::generate(void)
{
    for (int i = 0; i < mycontributions.size(); i++)
    {
        for (int j = 0; j < mycontributions[i].size(); j++)
            generate(i, j);
    }
}

void formulation::generatestiffnessmatrix(void)
{
    int i = 1;
    for (int j = 0; j < mycontributions[i].size(); j++)
        generate(i, j);
}

void formulation::generatedampingmatrix(void)
{
    int i = 2;
    for (int j = 0; j < mycontributions[i].size(); j++)
        generate(i, j);
}

void formulation::generatemassmatrix(void)
{
    int i = 3;
    for (int j = 0; j < mycontributions[i].size(); j++)
        generate(i, j);
}

void formulation::generaterhs(void)
{
    int i = 0;
    for (int j = 0; j < mycontributions[i].size(); j++)
        generate(i, j);
}


void formulation::generatein(int rhskcm, std::vector<int> contributionnumbers)
{
    for (int i = 0; i < contributionnumbers.size(); i++)
        generate(rhskcm, contributionnumbers[i]);
}

void formulation::generate(std::vector<int> contributionnumbers)
{
    for (int i = 0; i < mycontributions.size(); i++)
    {
        for (int j = 0; j < contributionnumbers.size(); j++)
            generate(i, contributionnumbers[j]);
    }
}

void formulation::generate(int contributionnumber)
{
    for (int i = 0; i < mycontributions.size(); i++)
            generate(i, contributionnumber);
}

densematrix formulation::getportrelationrhs(void)
{
    int expectednumrelations = mydofmanager->countports() - mydofmanager->countassociatedprimalports();

    densematrix rhsvals(expectednumrelations, 1, 0.0);
    double* vptr = rhsvals.getvalues();

    int actualnumrelations = 0;
    for (int i = 0; i < myportrelations.size(); i++)
    {
        if (myportrelations[i]->hasnoportterm())
            vptr[actualnumrelations] = -myportrelations[i]->evalnoportterm();
        
        actualnumrelations += myportrelations[i]->count();
    }

    if (expectednumrelations != actualnumrelations)
    {
        std::cout << "Error in 'formulation' object: expected " << expectednumrelations << " port relations to match the number of unknown ports provided but found " << actualnumrelations << std::endl;
        abort();
    }

    return rhsvals;
}

std::tuple<intdensematrix, intdensematrix, densematrix> formulation::getportrelations(int KCM)
{
    int expectednumrelations = mydofmanager->countports() - mydofmanager->countassociatedprimalports();
    
    std::vector< std::vector<int> > rinds(myportrelations.size());
    std::vector< std::vector<int> > cinds(myportrelations.size());
    std::vector< std::vector<double> > vals(myportrelations.size());

    int portnnz = 0;
    
    int actualnumrelations = 0;
    for (int i = 0; i < myportrelations.size(); i++)
    {
        std::vector<std::shared_ptr<rawport>> rps;
        myportrelations[i]->evalrelations(KCM, rps, rinds[i], vals[i]);
    
        int len = rps.size();
        cinds[i].resize(len);
    
        for (int j = 0; j < len; j++)
        {
            rinds[i][j] += actualnumrelations;
            cinds[i][j] = mydofmanager->getaddress(rps[j].get());
        }
    
        portnnz += len;
        actualnumrelations += myportrelations[i]->count();
    }
    
    if (expectednumrelations != actualnumrelations)
    {
        std::cout << "Error in 'formulation' object: expected " << expectednumrelations << " port relations to match the number of unknown ports provided but found " << actualnumrelations << std::endl;
        abort();
    }
    
    // Concatenate all:
    intdensematrix allrinds(portnnz, 1);
    intdensematrix allcinds(portnnz, 1);
    densematrix allvals(portnnz, 1);
    
    int* arptr = allrinds.getvalues();
    int* acptr = allcinds.getvalues();
    double* avptr = allvals.getvalues();
    
    int index = 0;
    for (int i = 0; i < myportrelations.size(); i++)
    {
        for (int j = 0; j < rinds[i].size(); j++)
        {
            arptr[index] = rinds[i][j];
            acptr[index] = cinds[i][j];
            avptr[index] = vals[i][j];
        
            index++;
        }
    }
    
    return std::make_tuple(allrinds, allcinds, allvals);
}


vec formulation::b(bool keepvector, bool dirichletandportupdate) { return rhs(keepvector, dirichletandportupdate); }
mat formulation::A(bool keepfragments) { return K(keepfragments); }

vec formulation::rhs(bool keepvector, bool dirichletandportupdate)
{
    if (myvec == NULL)
        myvec = std::shared_ptr<rawvec>(new rawvec(mydofmanager));
    
    vec output;   
    if (keepvector == false)
    {
        output = vec(myvec);
        myvec = NULL;
    }
    else
        output = vec(myvec).copy();
    
    if (dirichletandportupdate == true && isconstraintcomputation == false)
        output.updateconstraints(); 
        
    if (dirichletandportupdate == true)
    {
        densematrix portrhsvals = getportrelationrhs();
        output.setvalues(intdensematrix(portrhsvals.count(), 1, 0, 1), portrhsvals);
    }
    
    return output; 
}

mat formulation::K(bool keepfragments) { return getmatrix(0, keepfragments); }
mat formulation::C(bool keepfragments) { return getmatrix(1, keepfragments); }
mat formulation::M(bool keepfragments) { return getmatrix(2, keepfragments); }

mat formulation::getmatrix(int KCM, bool keepfragments, std::vector<intdensematrix> additionalconstraints)
{
    if (mymat[KCM] == NULL)
        mymat[KCM] = std::shared_ptr<rawmat>(new rawmat(mydofmanager));
        
    std::shared_ptr<rawmat> rawout = mymat[KCM]->extractaccumulated();
    
    if (keepfragments == false)
        mymat[KCM] = NULL;
        
    std::vector<bool> isconstr;
    if (isconstraintcomputation)
        isconstr = std::vector<bool>(mydofmanager->countdofs(), false);
    else
        isconstr = mydofmanager->isconstrained();
        
    for (int i = 0; i < additionalconstraints.size(); i++)
    {
        int* acptr = additionalconstraints[i].getvalues();
        for (int j = 0; j < additionalconstraints[i].count(); j++)
            isconstr[acptr[j]] = true;
    }

    if (KCM == 0)
    {
        std::pair<intdensematrix, intdensematrix> assocports = mydofmanager->findassociatedports();
        rawout->accumulate(assocports.first, assocports.second, densematrix(assocports.first.count(), 1, 1.0));
    }
    std::tuple<intdensematrix, intdensematrix, densematrix> portterms = getportrelations(KCM);
    rawout->accumulate(std::get<0>(portterms), std::get<1>(portterms), std::get<2>(portterms));
    
    rawout->process(isconstr); 
    rawout->clearfragments();
    
    return mat(rawout);
}

void formulation::solve(std::string soltype, bool diagscaling, std::vector<int> blockstoconsider)
{
    // Make sure the problem is of the form Ax = b:
    if (isdampingmatrixdefined() || ismassmatrixdefined())
    {
        std::cout << "Error in 'formulation' object: cannot solve with a damping/mass matrix (use a time resolution algorithm)" << std::endl;
        abort();  
    }
    
    // Remove leftovers (if any):
    myvec = NULL; mymat = {NULL, NULL, NULL};
    
    // Generate:
    if (blockstoconsider.size() == 1 && blockstoconsider[0] == -1)
        generate();
    else
        generate(blockstoconsider);
    // Solve:
    vec sol = sl::solve(this->A(), this->b(), soltype, diagscaling);

    // Save to fields:
    sl::setdata(sol);
}

densematrix Fgmultdirichlet(densematrix gprev)
{
    mat A = universe::ddmmats[0];
    formulation formul = universe::ddmformuls[0];
    
    vec rhs(formul);

    double* gprevptr = gprev.getvalues();

    std::shared_ptr<dtracker> dt = universe::mymesh->getdtracker();
    int numneighbours = dt->countneighbours();

    // Compute Ag:
    int pos = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int len = universe::ddmrecvinds[n].count();
        densematrix Bm = gprev.extractrows(pos, pos+len-1);
        rhs.setvalues(universe::ddmrecvinds[n], Bm);
        pos += len;
    }
        
    vec sol = sl::solve(A, rhs);

    // Create the artificial sources solution on the inner interface:
    std::vector<densematrix> Agmatssend(numneighbours), Agmatsrecv(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        Agmatssend[n] = sol.getvalues(universe::ddmsendinds[n]);
        Agmatsrecv[n] = densematrix(universe::ddmrecvinds[n].count(), 1);
    }
    sl::exchange(dt->getneighbours(), Agmatssend, Agmatsrecv);
    
    densematrix Ag(Agmatsrecv);
    
    // Calculate I - A*g:
    double* Agptr = Ag.getvalues();
    densematrix Fg(Ag.count(), 1);
    double* Fgptr = Fg.getvalues();
    
    for (int i = 0; i < Ag.count(); i++)
        Fgptr[i] = gprevptr[i] - Agptr[i];

    return Fg;
}

std::vector<double> formulation::allsolve(double relrestol, int maxnumit, std::string soltype, int verbosity)
{
    // Make sure the problem is of the form Ax = b:
    if (isdampingmatrixdefined() || ismassmatrixdefined())
    {
        std::cout << "Error in 'formulation' object: cannot solve with a damping/mass matrix (use a time resolution algorithm)" << std::endl;
        abort();  
    }
    
    wallclock clktot;

    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    if (numranks == 1)
    {
        solve(soltype, false);
        return {};
    }
    
    universe::ddmformuls = {*this};

    std::shared_ptr<dtracker> dt = universe::mymesh->getdtracker();
    int numneighbours = dt->countneighbours();

    // Get the rows from which to take the dofs to send as well as the rows at which to place the received dofs:
    std::vector<std::shared_ptr<rawfield>> rfs = mydofmanager->getfields();
    sl::mapdofs(mydofmanager, mydofmanager->getfields(), {true, true, true}, universe::ddmsendinds, universe::ddmrecvinds);
    
    std::vector<intdensematrix> interfaceinds = universe::ddmrecvinds;
    
    // Get all Dirichlet constraints set on the neighbours but not on this rank:
    std::vector<std::vector<intdensematrix>> dcdata = mydofmanager->discovernewconstraints(dt->getneighbours(), universe::ddmsendinds, universe::ddmrecvinds);

    // Unconstrained send and receive indexes:
    universe::ddmsendinds = dcdata[2];
    universe::ddmrecvinds = dcdata[3];
    
    // Get A and allow to reuse its factorization:
    generate();
    mat A = getmatrix(0, false, interfaceinds);
    A.reusefactorization();
    universe::ddmmats = {A};
    
    // Get the rhs of the physical sources contribution:
    vec bphysical = b();
    for (int n = 0; n < numneighbours; n++)
        bphysical.setvalues(dcdata[3][n], densematrix(dcdata[3][n].count(), 1, 0.0));
    
    // Set the value from the neighbour Dirichlet conditions:
    std::vector<densematrix> dirichletvalsforneighbours(numneighbours), dirichletvalsfromneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        dirichletvalsforneighbours[n] = bphysical.getvalues(dcdata[0][n]);
        dirichletvalsfromneighbours[n] = densematrix(dcdata[1][n].count(), 1);
    }
    sl::exchange(dt->getneighbours(), dirichletvalsforneighbours, dirichletvalsfromneighbours);
    for (int n = 0; n < numneighbours; n++)
        bphysical.setvalues(dcdata[1][n], dirichletvalsfromneighbours[n]);

    // Get the physical sources solution:
    vec w = sl::solve(A, bphysical, soltype);
    
    // Create the gmres rhs 'B' from the physical sources solution on the inner interface:
    std::vector<densematrix> Bmatssend(numneighbours), Bmatsrecv(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        Bmatssend[n] = w.getvalues(universe::ddmsendinds[n]);
        Bmatsrecv[n] = densematrix(universe::ddmrecvinds[n].count(), 1);
    }
    sl::exchange(dt->getneighbours(), Bmatssend, Bmatsrecv);
    
    densematrix B(Bmatsrecv);
    
    // Initial value of the artificial sources solution:
    densematrix vi(B.countrows(), B.countcolumns(), 0.0);
    
    // Gmres iteration:
    std::vector<double> resvec = sl::gmres(Fgmultdirichlet, B, vi, relrestol, maxnumit, verbosity*(rank == 0));
    int numits = resvec.size()-1;
    if (verbosity > 0 && rank == 0)
    {
        if (numits < maxnumit)
            std::cout << "gmres converged after " << numits << " iterations" << std::endl;
        else
            std::cout << "gmres could not converge to the requested " << relrestol << " relative tolerance (could only reach " << resvec[numits] << ")" << std::endl;
    }

    // Compute the total solution (physical + artificial):
    int pos = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int len = Bmatsrecv[n].count();
        densematrix Bm = vi.extractrows(pos, pos+len-1);
        bphysical.setvalues(universe::ddmrecvinds[n], Bm);
        pos += len;
    }

    vec totalsol = sl::solve(A, bphysical);
    sl::setdata(totalsol);
    
    if (verbosity > 0)
    {
        long long int allndfs = allcountdofs();
        if (rank == 0)
            clktot.print("DDM solve for "+std::to_string(allndfs)+" dofs took");
    }
    
    universe::clearddmcontainers();
    
    return resvec;
}

densematrix Fgmultrobin(densematrix gprev)
{
    mat A = universe::ddmmats[0];
    formulation formul = universe::ddmformuls[0];
    std::vector<std::vector<int>> artificialterms = universe::ddmints;

    vec rhs(formul);

    double* gprevptr = gprev.getvalues();

    std::shared_ptr<dtracker> dt = universe::mymesh->getdtracker();
    int numneighbours = dt->countneighbours();

    // Compute Ag:
    int pos = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int len = universe::ddmrecvinds[n].count();
        densematrix Bm = gprev.extractrows(pos, pos+len-1);
        rhs.setvalues(universe::ddmrecvinds[n], Bm, "add");
        pos += len;
    }
        
    vec sol = sl::solve(A, rhs);
    sl::setdata(sol);

    // Create the artificial sources solution on the inner interface:
    std::vector<densematrix> Agmatssend(numneighbours), Agmatsrecv(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        formul.generatein(0, artificialterms[n]);
        vec gartificial = formul.b(false, false);
    
        Agmatssend[n] = gartificial.getvalues(universe::ddmsendinds[n]);
        Agmatsrecv[n] = densematrix(universe::ddmrecvinds[n].count(), 1);
    }
    sl::exchange(dt->getneighbours(), Agmatssend, Agmatsrecv);
    
    densematrix Ag(Agmatsrecv);
    
    // Calculate I - A*g:
    double* Agptr = Ag.getvalues();
    densematrix Fg(Ag.count(), 1);
    double* Fgptr = Fg.getvalues();
    
    for (int i = 0; i < Ag.count(); i++)
        Fgptr[i] = gprevptr[i] - Agptr[i];

    return Fg;
}

std::vector<double> formulation::allsolve(std::vector<int> formulterms, std::vector<std::vector<int>> physicalterms, std::vector<std::vector<int>> artificialterms, double relrestol, int maxnumit, std::string soltype, int verbosity)
{
    // Make sure the problem is of the form Ax = b:
    if (isdampingmatrixdefined() || ismassmatrixdefined())
    {
        std::cout << "Error in 'formulation' object: cannot solve with a damping/mass matrix (use a time resolution algorithm)" << std::endl;
        abort();  
    }
    
    wallclock clktot;

    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    if (numranks == 1)
    {
        solve(soltype, false, formulterms);
        return {};
    }
    
    universe::ddmints = artificialterms;
    universe::ddmformuls = {*this};

    std::shared_ptr<dtracker> dt = universe::mymesh->getdtracker();
    int numneighbours = dt->countneighbours();

    // Get the rows from which to take the dofs to send as well as the rows at which to place the received dofs:
    std::vector<std::shared_ptr<rawfield>> rfs = mydofmanager->getfields();
    sl::mapdofs(mydofmanager, mydofmanager->getfields(), {true, true, true}, universe::ddmsendinds, universe::ddmrecvinds);
    
    // Get all Dirichlet constraints set on the neighbours but not on this rank:
    std::vector<std::vector<intdensematrix>> dcdata = mydofmanager->discovernewconstraints(dt->getneighbours(), universe::ddmsendinds, universe::ddmrecvinds);

    // Unconstrained send and receive indexes:
    universe::ddmsendinds = dcdata[2];
    universe::ddmrecvinds = dcdata[3];
    
    // Get A and allow to reuse its factorization:
    generate(formulterms);
    for (int n = 0; n < numneighbours; n++)
        generatein(1, physicalterms[n]); // S term in A
    mat A = getmatrix(0, false, {intdensematrix(dcdata[1])});
    A.reusefactorization();
    universe::ddmmats = {A};
    
    // Get the rhs of the physical sources contribution:
    vec bphysical = b();
    
    // Set the value from the neighbour Dirichlet conditions:
    std::vector<densematrix> dirichletvalsforneighbours(numneighbours), dirichletvalsfromneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        dirichletvalsforneighbours[n] = bphysical.getvalues(dcdata[0][n]);
        dirichletvalsfromneighbours[n] = densematrix(dcdata[1][n].count(), 1);
    }
    sl::exchange(dt->getneighbours(), dirichletvalsforneighbours, dirichletvalsfromneighbours);
    for (int n = 0; n < numneighbours; n++)
        bphysical.setvalues(dcdata[1][n], dirichletvalsfromneighbours[n]);

    // Get the physical sources solution:
    vec w = sl::solve(A, bphysical, soltype);
    sl::setdata(w);
    
    // Create the gmres rhs 'B' from the physical sources solution on the inner interface:
    std::vector<densematrix> Bmatssend(numneighbours), Bmatsrecv(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        generatein(0, physicalterms[n]);
        vec gphysical = b(false, false);
    
        Bmatssend[n] = gphysical.getvalues(universe::ddmsendinds[n]);
        Bmatsrecv[n] = densematrix(universe::ddmrecvinds[n].count(), 1);
    }
    sl::exchange(dt->getneighbours(), Bmatssend, Bmatsrecv);
    
    densematrix B(Bmatsrecv);
    
    // Initial value of the artificial sources solution:
    densematrix vi(B.countrows(), B.countcolumns(), 0.0);
    
    // Gmres iteration:
    std::vector<double> resvec = sl::gmres(Fgmultrobin, B, vi, relrestol, maxnumit, verbosity*(rank == 0));
    int numits = resvec.size()-1;
    if (verbosity > 0 && rank == 0)
    {
        if (numits < maxnumit)
            std::cout << "gmres converged after " << numits << " iterations" << std::endl;
        else
            std::cout << "gmres could not converge to the requested " << relrestol << " relative tolerance (could only reach " << resvec[numits] << ")" << std::endl;
    }

    // Compute the total solution (physical + artificial):
    int pos = 0;
    for (int n = 0; n < numneighbours; n++)
    {
        int len = Bmatsrecv[n].count();
        densematrix Bm = vi.extractrows(pos, pos+len-1);
        bphysical.setvalues(universe::ddmrecvinds[n], Bm, "add");
        pos += len;
    }

    vec totalsol = sl::solve(A, bphysical);
    sl::setdata(totalsol);
    
    if (verbosity > 0)
    {
        long long int allndfs = allcountdofs();
        if (rank == 0)
            clktot.print("DDM solve for "+std::to_string(allndfs)+" dofs took");
    }
    
    universe::clearddmcontainers();
    
    return resvec;
}

