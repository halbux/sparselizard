#include "formulation.h"


formulation::formulation(void) { mydofmanager = shared_ptr<dofmanager>(new dofmanager); }

void formulation::operator+=(integration integrationobject)
{
    int integrationphysreg = integrationobject.getphysicalregion();
    int elementdimension = universe::mymesh->getphysicalregions()->get(integrationphysreg)->getelementdimension();
    int integrationorderdelta = integrationobject.getintegrationorderdelta();
    expression myexpression = integrationobject.getexpression();

    myexpression.expand();
    
    // The element dimension is required to decompose space derivatives 
    // into the ki, eta and phi derivatives in the reference element.
    std::vector< std::vector<std::vector<shared_ptr<operation>>> > coeffdoftf = myexpression.extractdoftfpolynomial(elementdimension);

    std::vector<std::vector<shared_ptr<operation>>> coeffs = coeffdoftf[0]; 
    std::vector<std::vector<shared_ptr<operation>>> dofs = coeffdoftf[1];
    std::vector<std::vector<shared_ptr<operation>>> tfs = coeffdoftf[2];

    // Loop on all slices:
    for (int slice = 0; slice < tfs.size(); slice++)
    {
        // In a given slice all dof and tf fields are the same, have the 
        // same applied time derivatives and are selected on a same 
        // physical region for all entries.
        shared_ptr<rawfield> doffield = dofs[slice][0]->getfieldpointer();
        shared_ptr<rawfield> tffield = tfs[slice][0]->getfieldpointer();
        
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
        if (doffield != NULL)
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

void formulation::generate(int m, int contributionnumber)
{
    if (contributionnumber < 0)
    {
        std::cout << "Error in 'formulation' object: cannot generate a negative contribution number" << std::endl;
        abort();
    }
    
    // Make sure the contribution number exists:
    if (contributionnumber >= mycontributions[m].size() || mycontributions[m][contributionnumber].size() == 0)
        return;
        
    if (m == 0 && myvec == NULL)
        myvec = shared_ptr<rawvec>(new rawvec(mydofmanager));
    if (m > 0 && mymat[m-1] == NULL)
        mymat[m-1] = shared_ptr<rawmat>(new rawmat(mydofmanager));

    std::vector<contribution> contributionstogenerate = mycontributions[m][contributionnumber];
	for (int i = 0; i < contributionstogenerate.size(); i++)
    {
        if (m == 0)
            contributionstogenerate[i].generate(myvec, NULL, not(isconstraintcomputation));
        else
            contributionstogenerate[i].generate(NULL, mymat[m-1], not(isconstraintcomputation));
    }
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


vec formulation::b(bool keepvector) { return rhs(keepvector); }
mat formulation::A(bool keepfragments) { return K(keepfragments); }

vec formulation::rhs(bool keepvector)
{
    if (myvec == NULL)
        myvec = shared_ptr<rawvec>(new rawvec(mydofmanager));
    
    vec output;   
    if (keepvector == false)
    {
        output = vec(myvec);
        myvec = NULL;
    }
    else
        output = vec(myvec).copy();
    
	// Set the gauged indexes to zero:
	intdensematrix gaugedindexes = mydofmanager->getgaugedindexes();
	int numgaugedindexes = gaugedindexes.countrows()*gaugedindexes.countcolumns();
	if (numgaugedindexes > 0)
		output.getpointer()->setvalues(gaugedindexes, densematrix(gaugedindexes.countrows(),gaugedindexes.countcolumns(), 0.0));

    if (isconstraintcomputation == false)
        output.updateconstraints();  
    
    return output; 
}

mat formulation::K(bool keepfragments) { return getmatrix(0, keepfragments); }
mat formulation::C(bool keepfragments) { return getmatrix(1, keepfragments); }
mat formulation::M(bool keepfragments) { return getmatrix(2, keepfragments); }

mat formulation::getmatrix(int KCM, bool keepfragments)
{
    if (mymat[KCM] == NULL)
        mymat[KCM] = shared_ptr<rawmat>(new rawmat(mydofmanager));
        
    shared_ptr<rawmat> rawout = mymat[KCM]->extractaccumulated();
    
    if (keepfragments == false)
    	mymat[KCM]->clearfragments();
    	
        
	// Set the gauged indexes to all zero:
	rawout->gauge();

    // Add the constraint diagonal ones to the matrix (if any):
    if (isconstraintcomputation == false)
    {
        int numconstraineddofs = mydofmanager->countconstraineddofs();
        if (numconstraineddofs > 0)
        {
            intdensematrix constrainedindexes = mydofmanager->getconstrainedindexes();
            // Create a vector of same length full of double ones:
            densematrix ones(1, numconstraineddofs, 1);

            rawout->accumulate(constrainedindexes, constrainedindexes, ones);  
        }
    }
    // Add the gauge diagonal ones to the matrix (if any):
    int numgaugeddofs = mydofmanager->countgaugeddofs();
    if (numgaugeddofs > 0)
    {
        intdensematrix gaugedindexes = mydofmanager->getgaugedindexes();
        // Create a vector of same length full of double ones:
        densematrix ones(1, numgaugeddofs, 1);

        rawout->accumulate(gaugedindexes, gaugedindexes, ones);  
    }

    rawout->process(); 
    rawout->clearfragments();
    
    return mat(rawout);
}


