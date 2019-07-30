#include "rawfield.h"


rawfield::rawfield(std::string fieldtypename, const std::vector<int> harmonicnumbers, bool ismultiharm)
{
    multiharmonic = ismultiharm;
    
    // Treat the coordinate fields:
    if (fieldtypename == "x" || fieldtypename == "y" || fieldtypename == "z")
    {            
        myname = fieldtypename;
        mytypename = fieldtypename;
        // A coordinate field can only be constant (i.e. on harmonic 1). 
        if (harmonicnumbers.size() != 1 || harmonicnumbers[0] != 1)
        {
            std::cout << "Error in 'rawfield' object: rawfield type " << fieldtypename << " cannot be multiharmonic (must be constant)" << std::endl;
            abort();
        }
        return;
    }

    // Make sure the mesh has been loaded before defining non-coordinate fields:
    if (universe::mymesh == NULL)
    {
        std::cout << "Error in 'rawfield' object: first load mesh before defining a field that is not the x, y or z coordinate" << std::endl;
        abort();
    }
    
    // Edge shape functions cannot be used with axisymmetry for now:
	if (universe::isaxisymmetric && fieldtypename == "hcurl")
	{
		std::cout << "Error in 'rawfield' object: using hcurl (edge) shape functions with axisymmetry is not allowed." << std::endl;
		std::cout << "Consider using an alternative formulation with h1 (nodal) shape functions." << std::endl;
		abort();
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

    // Make sure the subfield type is scalar:
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(0, mytypename);
    if (numberofsubfields > 1 && myformfunction->countcomponents() != 1)
    {
        std::cout << "Error in 'rawfield' object: " << fieldtypename << " is not a valid field type. Only scalar form functions can have a trailing 'xy' or 'xyz'" << std::endl;
        abort();
    }

    // For fields with subfields:
    if (numberofsubfields > 1)
    {
        for (int i = 0; i < numberofsubfields; i++)
            mysubfields.push_back({ shared_ptr<rawfield>(new rawfield(mytypename, harmonicnumbers, ismultiharm)) });
    }
    else
    {
        if (harmonicnumbers.size() != 0)
        {
            myharmonics.resize(*max_element(harmonicnumbers.begin(), harmonicnumbers.end())+1);
            for (int h = 0; h < harmonicnumbers.size(); h++)
                myharmonics[harmonicnumbers[h]] = { shared_ptr<rawfield>(new rawfield(mytypename, {}, ismultiharm)) };
        }    
        else
        {
            // Set an order 1 interpolation order by default:
            interpolationorder = std::vector<int>( (universe::mymesh->getdisjointregions())->count(), 1);
            // Set all unconstrained by default:
            myconstraints = std::vector<shared_ptr<integration>>( (universe::mymesh->getdisjointregions())->count(), NULL);
            
            myconditionalconstraints = std::vector<std::vector<expression>>( (universe::mymesh->getdisjointregions())->count(), std::vector<expression>(0));

            isitgauged = std::vector<bool>( (universe::mymesh->getdisjointregions())->count(), false);

            mycoefmanager = std::shared_ptr<coefmanager>(new coefmanager(mytypename));
        }
    }
    return;
}

int rawfield::countcomponents(void)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
        return 1;
        
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(0, mytypename);
    return std::max(myformfunction->countcomponents(), (int) mysubfields.size());
}

int rawfield::countformfunctioncomponents(void)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
        return 1;
    
    // Get the number of form function components:
    std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(0, mytypename);
    return myformfunction->countcomponents();
}

std::vector<int> rawfield::getharmonics(void)
{
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
    if (mysubfields.size() == 0)
    {
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() == 1)
                return h;
        }
    }
    else
        return mysubfields[0][0]->getfirstharmonic();
}

bool rawfield::isharmonicincluded(int harmonic)
{
    if (mysubfields.size() == 0)
        return (harmonic == 1 && myharmonics.size() == 0) || (harmonic > 0 && harmonic < myharmonics.size() && myharmonics[harmonic].size() > 0);
    else
        return mysubfields[0][0]->isharmonicincluded(harmonic);
}

void rawfield::printharmonics(void)
{
    if (multiharmonic == false)
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

void rawfield::print(void)
{
    if (myname.size() == 0)
        std::cout << "field";
    else
        std::cout << myname;
}

void rawfield::setorder(int physreg, int interpolorder)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z" || mytypename == "one")
    {
        std::cout << "Error in 'rawfield' object: cannot choose the interpolation order for the x, y, z coordinate or for 'one' type fields" << std::endl;
        abort();
    }
        
    // Set the interpolation order on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setorder(physreg, interpolorder);
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setorder(physreg, interpolorder);
    }
        
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
    {
        // Consider ALL disjoint regions in the physical region with (-1):
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);
        
        for (int i = 0; i < selecteddisjregs.size(); i++)
        {
            interpolationorder[selecteddisjregs[i]] = interpolorder;
            mycoefmanager->fitinterpolationorder(selecteddisjregs[i], interpolorder);
        }
    }
}

void rawfield::setvalue(int physreg, int numfftharms, expression* meshdeform, expression input, int extraintegrationdegree)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        std::cout << "Error in 'rawfield' object: cannot set the value for the x, y or z coordinate" << std::endl;
        abort();
    }
	if (input.countcolumns() != 1 || input.countrows() != countcomponents())
    {
        std::cout << "Error in 'rawfield' object: the rawfield value must be set with a " << countcomponents() << "x1 expression" << std::endl;
        abort();
    }
    
    // Set the values on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setvalue(physreg, numfftharms, meshdeform, input.at(i,0), extraintegrationdegree);

    if (mysubfields.size() == 0)
    {
    	field thisfield(getpointer());
    
    	// Compute the projection of the expression (skip for a zero expression):
    	formulation projectedvalue;
    	if (meshdeform == NULL)
    		projectedvalue += integration(physreg, numfftharms, mathop::transpose(mathop::dof(thisfield))*mathop::tf(thisfield) - mathop::transpose(mathop::tf(thisfield))*input, extraintegrationdegree);
    	else
    		projectedvalue += integration(physreg, numfftharms, *meshdeform, mathop::transpose(mathop::dof(thisfield))*mathop::tf(thisfield) - mathop::transpose(mathop::tf(thisfield))*input, extraintegrationdegree);
    	// Define an all-zero vector:
    	vec solvec(projectedvalue);
		if (input.iszero() == false)
		{
			projectedvalue.generate();
			solvec = mathop::solve(projectedvalue.A(), projectedvalue.b());
		}
    	
    	thisfield.setdata(physreg, solvec);
    }
}

void rawfield::setvalue(int physreg)
{
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

void rawfield::setconstraint(int physreg, expression input, int extraintegrationdegree)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        std::cout << "Error in 'rawfield' object: cannot constrain the x, y or z coordinate" << std::endl;
        abort();
    }
    if (input.countcolumns() != 1 || input.countrows() != countcomponents())
    {
        std::cout << "Error in 'rawfield' object: the rawfield must be constrained using a " << countcomponents() << "x1 expression" << std::endl;
        abort();
    }
        
    // Set the constraints on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setconstraint(physreg, input.at(i,0), extraintegrationdegree);
    // Set the constraints on the harmonics:
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setconstraint(physreg, input, extraintegrationdegree);
    }
    
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
    {
        // Consider ALL disjoint regions in the physical region with (-1):
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);

        field thisfield(getpointer());
        std::shared_ptr<integration> constraintcomputation(new integration(physreg, mathop::transpose(mathop::dof(thisfield))*mathop::tf(thisfield) - mathop::transpose(mathop::tf(thisfield))*input, extraintegrationdegree));
        
        if (input.iszero())
            constraintcomputation->isprojectionofzero = true;
        
        for (int i = 0; i < selecteddisjregs.size(); i++)
            myconstraints[selecteddisjregs[i]] = constraintcomputation;
    }
}

void rawfield::setconstraint(int physreg)
{
    switch (countcomponents())
    {
        case 1:
            setconstraint(physreg, 0);
            break;
        case 2:
            setconstraint(physreg, expression(2,1,{0,0}));
            break;
        case 3:
            setconstraint(physreg, expression(3,1,{0,0,0}));
            break;
    }
}

void rawfield::setconditionalconstraint(int physreg, expression condexpr, expression valexpr)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        std::cout << "Error in 'rawfield' object: cannot constrain the x, y or z coordinate" << std::endl;
        abort();
    }
    if (valexpr.countcolumns() != 1 || valexpr.countrows() != countcomponents())
    {
        std::cout << "Error in 'rawfield' object: the rawfield must be constrained using a " << countcomponents() << "x1 expression" << std::endl;
        abort();
    }
    if (condexpr.countcolumns() != 1 || condexpr.countrows() != 1)
    {
        std::cout << "Error in 'rawfield' object: expected a scalar condition for the conditional constraint" << std::endl;
        abort();
    }
        
    // Set the constraints on the sub fields:
    for (int i = 0; i < mysubfields.size(); i++)
        mysubfields[i][0]->setconditionalconstraint(physreg, condexpr, valexpr.at(i,0));
    // Set the constraints on the harmonics:
    for (int i = 0; i < myharmonics.size(); i++)
    {
        if (myharmonics[i].size() > 0)
            myharmonics[i][0]->setconditionalconstraint(physreg, condexpr, valexpr);
    }
    
    if (mysubfields.size() == 0 && myharmonics.size() == 0)
    {
        // Consider only the NODAL disjoint regions in the physical region with (0):
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(0);

        for (int i = 0; i < selecteddisjregs.size(); i++)
            myconditionalconstraints[selecteddisjregs[i]] = {condexpr, valexpr};
    }
}

void rawfield::setgauge(int physreg)
{
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
        std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);

		for (int i = 0; i < selecteddisjregs.size(); i++)
			isitgauged[selecteddisjregs[i]] = true;
    }
}

void rawfield::setspanningtree(spanningtree spantree)
{
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
    	myspanningtree = {spantree};
}

spanningtree* rawfield::getspanningtree(void)
{
	if (myspanningtree.size() == 1)
		return &myspanningtree[0];
	else
	{
		std::cout << "Error in 'rawfield' object: spanning tree was not provided to rawfield" << std::endl;
		abort();
	}
}

void rawfield::setdata(int physreg, vectorfieldselect myvec, std::string op)
{
    // Extract the info from the vector with selected field:
    shared_ptr<rawfield> selectedrawfield = myvec.getrawfield();
    shared_ptr<rawvec> selectedvec = myvec.getrawvector();
        
    // The raw fields must be of the same type:
    if (mytypename != selectedrawfield->mytypename)
    {
        std::cout << "Error in 'rawfield' object: .setdata can only transfer data between fields of same type" << std::endl;
        abort();
    }
    // The raw fields must have a same number of subfields:
    if (mysubfields.size() != selectedrawfield->mysubfields.size())
    {
        std::cout << "Error in 'rawfield' object: .setdata can only transfer data from fields with same number of subfields" << std::endl;
        abort();
    }
    
    // Get the data of every subfield:
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
            mysubfields[i][0]->setdata(physreg, vectorfieldselect(selectedvec, selectedrawfield->mysubfields[i][0]), op);
        return;
    }
    
    // The raw fields must include the same harmonic numbers:
    if (getharmonics() != selectedrawfield->getharmonics())
    {
        std::cout << "Error in 'rawfield' object: .setdata can only transfer data from fields with same harmonic numbers" << std::endl;
        abort();
    }
    
    // Get the data of every harmonic:
    if (myharmonics.size() > 0)
    {
        for (int h = 0; h < myharmonics.size(); h++)
        {
            if (myharmonics[h].size() > 0)
                myharmonics[h][0]->setdata(physreg, vectorfieldselect(selectedvec, selectedrawfield->harmonic(h)), op);
        }
        return;
    }
    // Extract the actual field from non-multiharmonic fields:
    if (selectedrawfield->multiharmonic == false)
        selectedrawfield = selectedrawfield->harmonic(1);
    
    // Get the data for a single field.
    
    // Get ALL disjoint regions in the physical region:
    std::vector<int> selecteddisjregs;
    if (physreg != -1)
        selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);
    else
    {
        std::shared_ptr<dofmanager> dofmngr = selectedvec->getdofmanager();
        dofmngr->selectfield(shared_from_this());
        selecteddisjregs = dofmngr->getdisjointregionsofselectedfield();
    }

    for (int i = 0; i < selecteddisjregs.size(); i++)
    {
        int disjreg = selecteddisjregs[i];

        // In case the order of this raw field is higher than the order of 
        // the selected raw field we have to set to zero the higher orders.
        if (op == "set" && interpolationorder[disjreg] > selectedrawfield->interpolationorder[disjreg])
        {
            // Decrease the order to forget the higher orders...
            mycoefmanager->fitinterpolationorder(disjreg, selectedrawfield->interpolationorder[disjreg]);
            // ... then reset the previous order:
            mycoefmanager->fitinterpolationorder(disjreg, interpolationorder[disjreg]);
        }

        int elementtypenumber = (universe::mymesh->getdisjointregions())->getelementtypenumber(disjreg);
        int elementdimension = (universe::mymesh->getdisjointregions())->getelementdimension(disjreg);
        int numelem = (universe::mymesh->getdisjointregions())->countelements(disjreg);

        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, mytypename);
        // The interpolation order for this field and the selected fields might be different.
        int numformfunctionsperelement = std::min(myformfunction->count(interpolationorder[disjreg], elementdimension, 0), myformfunction->count(selectedrawfield->interpolationorder[disjreg], elementdimension, 0));

        for (int ff = 0; ff < numformfunctionsperelement; ff++)
        {
            densematrix values = selectedvec->getvalues(selectedrawfield, disjreg, ff);
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

void rawfield::transferdata(int physreg, vectorfieldselect myvec, std::string op)
{
    // Extract the info from the vector with selected field:
    shared_ptr<rawfield> selectedrawfield = myvec.getrawfield();
    shared_ptr<rawvec> selectedvec = myvec.getrawvector();
        
    // The raw fields must be of the same type:
    if (mytypename != selectedrawfield->mytypename)
    {
        std::cout << "Error in 'rawfield' object: .transferdata can only transfer data between fields of same type" << std::endl;
        abort();
    }
    // The raw fields must have a same number of subfields:
    if (mysubfields.size() != selectedrawfield->mysubfields.size())
    {
        std::cout << "Error in 'rawfield' object: .transferdata can only transfer data from fields with same number of subfields" << std::endl;
        abort();
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
        std::cout << "Error in 'rawfield' object: .transferdata can only transfer data from fields with same harmonic numbers" << std::endl;
        abort();
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
    if (selectedrawfield->multiharmonic == false)
        selectedrawfield = selectedrawfield->harmonic(1);
    
    // Transfer the data for a single field.
    
    // Get ALL disjoint regions in the physical region:
    std::vector<int> selecteddisjregs;
    if (physreg != -1)
        selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);
    else
    {
        std::shared_ptr<dofmanager> dofmngr = selectedvec->getdofmanager();
        dofmngr->selectfield(shared_from_this());
        selecteddisjregs = dofmngr->getdisjointregionsofselectedfield();
    }

    for (int i = 0; i < selecteddisjregs.size(); i++)
    {
        int disjreg = selecteddisjregs[i];

        int elementtypenumber = (universe::mymesh->getdisjointregions())->getelementtypenumber(disjreg);
        int elementdimension = (universe::mymesh->getdisjointregions())->getelementdimension(disjreg);
        int numelem = (universe::mymesh->getdisjointregions())->countelements(disjreg);

        std::shared_ptr<hierarchicalformfunction> myformfunction = selector::select(elementtypenumber, mytypename);
        // The interpolation order for this field and the selected fields might be different.
        int numformfunctionsinoriginfield = myformfunction->count(interpolationorder[disjreg], elementdimension, 0);
        int numformfunctionsperelement = myformfunction->count(selectedrawfield->interpolationorder[disjreg], elementdimension, 0);

        for (int ff = 0; ff < numformfunctionsperelement; ff++)
        {
            densematrix values(1,numelem,0.0);
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

shared_ptr<rawfield> rawfield::comp(int component)
{   
    // If there is a single component and the first one is requested:
    if (countcomponents() == 1 && component == 0)
        return shared_from_this();
    
    if (countformfunctioncomponents() > 1)
    {
        std::cout << "Error in 'rawfield' object: cannot get a component for vector fields with no subfields (e.g. hcurl)" << std::endl;
        abort();
    }
    if (component > mysubfields.size())
    {
        std::cout << "Error in 'rawfield' object: cannot get component " << component << " from a " << mysubfields.size() << " components field" << std::endl;
        abort();
    }
        
    return mysubfields[component][0];
}

shared_ptr<rawfield> rawfield::harmonic(const std::vector<int> harmonicnumbers)
{
    if (mysubfields.size() != 0)
    {
        shared_ptr<rawfield> harmsrawfield(new rawfield());
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
                std::cout << "Error in 'rawfield' object: in .harmonic cannot get harmonic " << harmonicnumbers[0] << " (does not exist)" << std::endl; 
                abort();
            }
        }
        
        // In case we want several harmonics a new raw field container is created.
        int maxharmnum = *max_element(harmonicnumbers.begin(), harmonicnumbers.end());

        shared_ptr<rawfield> harmsrawfield(new rawfield());
        *harmsrawfield = *this;

        // Set a brand new harmonic vector:
        harmsrawfield->myharmonics = std::vector<std::vector<shared_ptr<rawfield>>>(maxharmnum+1, std::vector<shared_ptr<rawfield>>(0));
        
        // Add only the requested harmonics:
        for (int i = 0; i < harmonicnumbers.size(); i++)
        {
            if (harmonicnumbers[i] < myharmonics.size() && myharmonics[harmonicnumbers[i]].size() > 0)
                harmsrawfield->myharmonics[harmonicnumbers[i]] = {myharmonics[harmonicnumbers[i]][0]};
            else
            {
                std::cout << "Error in 'rawfield' object: in .harmonic cannot get harmonic " << harmonicnumbers[i] << " (does not exist)" << std::endl; 
                abort();
            }
        }
        return harmsrawfield;
    }
        
    // In case there is no harmonic it is a constant (harmonic 1).
    if (harmonicnumbers.size() == 1 && harmonicnumbers[0] == 1)
        return shared_from_this();
    else
    {
        std::cout << "Error in 'rawfield' object: in .harmonic cannot get harmonic in constant field (does not exist)" << std::endl; 
        abort();
    }
}

bool rawfield::isgauged(int disjreg) 
{ 
	return isitgauged[disjreg];
}

int rawfield::getinterpolationorder(int disjreg) 
{ 
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
        return 1;
        
    if (mysubfields.size() == 0)
    {
        if (myharmonics.size() == 0)
            return interpolationorder[disjreg];
        else
        {
            errornotsameinterpolationorder(disjreg);
            return myharmonics[getfirstharmonic()][0]->interpolationorder[disjreg];
        }
    }
    else
    { 
        std::cout << "Error in 'rawfield' object: cannot get the interpolation order of a field with subfields" << std::endl;
        abort();
    }
}

void rawfield::errornotsameinterpolationorder(int disjreg)
{
    if (mysubfields.size() == 0)
    {
        if (myharmonics.size() > 0)
        {
            std::vector<int> harms = getharmonics();
            for (int h = 0; h < harms.size(); h++)
            {
                if (myharmonics[harms[h]][0]->interpolationorder[disjreg] == myharmonics[harms[0]][0]->interpolationorder[disjreg])
                    continue;
                else
                {
                    std::cout << "Error in 'rawfield' object: the interpolation order must be the same for all harmonics" << std::endl;
                    abort();
                }
            }
        }
    }
    else
    { 
        std::cout << "Error in 'rawfield' object: cannot call 'errornotsameinterpolationorder' on a field with subfields" << std::endl;
        abort();
    }
}


std::vector<std::pair<std::vector<int>, shared_ptr<rawfield>>> rawfield::getallsons(void)
{
    std::vector<std::pair<std::vector<int>, shared_ptr<rawfield>>> output = {};
    
    if (mysubfields.size() > 0)
    {
        for (int i = 0; i < mysubfields.size(); i++)
        {
            if (mysubfields[i].size() > 0)
            {
                std::vector<std::pair<std::vector<int>, shared_ptr<rawfield>>> cur = mysubfields[i][0]->getallsons();
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
                std::vector<std::pair<std::vector<int>, shared_ptr<rawfield>>> cur = myharmonics[h][0]->getallsons();
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

void rawfield::writeraw(int physreg, std::string filename, bool isbinary)
{
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        std::cout << "Error in 'rawfield' object: cannot write field type '" << mytypename << "' to raw format" << std::endl;
        abort();
    }
    
    int numsubfields = countsubfields();
    std::vector<int> harmoniclist = getharmonics();
    int numharms = harmoniclist.size();

    hierarchicalformfunction myhff;
    int fieldtypenum = myhff.gettypenumber(mytypename);

    // Get all disjoint regions in the physical region:
    std::vector<int> selecteddisjregs = ((universe::mymesh->getphysicalregions())->get(physreg))->getdisjointregions(-1);
    int numdisjregs = selecteddisjregs.size();
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();

    // Get all sons:
    std::vector<std::pair<std::vector<int>, shared_ptr<rawfield>>> allsons = getallsons();
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


    int doubledatasize = 0;
    for (int i = 0; i < numsons; i++)
    {
        shared_ptr<rawfield> curson = allsons[i].second;
        shared_ptr<coefmanager> curcoefmanager = curson->mycoefmanager;
    
        // For every son provide info for each disjoint region as {order, doublevecrangebegin}:
        for (int d = 0; d < numdisjregs; d++)
        {
            int disjreg = selecteddisjregs[d];
        
            int interpolorder = curson->interpolationorder[disjreg];
         
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

void rawfield::loadraw(std::string filename, bool isbinary)
{
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();
    
    
    ///// Load the data from disk:
    std::vector<int> intdata;
    std::vector<double> doubledata;

    iointerface::load(filename, intdata, doubledata, isbinary);
    
    
    ///// Extract the file format number:
    int fileformatnum = intdata[0];
    
    
    ///// Extract the field type number and make sure the field type names match:
    int fieldtypenum = intdata[1];
    hierarchicalformfunction myhff;
    std::string fieldtypename = myhff.gettypename(fieldtypenum);

    if (mytypename != fieldtypename)
    {
        std::cout << "Error in 'rawfield' object: trying to load a '" << fieldtypename << "' type field from file '" << filename << "' to a '" << mytypename << "' type" << std::endl;
        abort();
    }
    
    
    ///// Extract the number of subfields and make sure they match:
    int numsubfields = intdata[2];
    if (numsubfields != countsubfields())
    {
        std::cout << "Error in 'rawfield' object: trying to load a field with " << numsubfields << " subfields from file '" << filename << "' to one with " << countsubfields() << " subfields" << std::endl;
        abort();
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
        std::cout << "Error in 'rawfield' object: harmonic list in field loaded from file '" << filename << "' does not match current field." << std::endl << std::endl;

        std::cout << "Harmonics in loaded field (" << harmoniclist.size() << "): ";
        for (int i = 0; i < harmoniclist.size(); i++)
            std::cout << harmoniclist[i] << " ";
        std::cout << std::endl;
        
        std::cout << "Harmonics in current field (" << harmsinthisfield.size() << "): ";
        for (int i = 0; i < harmsinthisfield.size(); i++)
            std::cout << harmsinthisfield[i] << " ";
        std::cout << std::endl << std::endl;
        abort();
    }
    
    
    ///// Extract the list of disjoint regions and the number of elements in them from intdata, make sure the meshes match:
    int numdisjregs = intdata[indexinintvec];
    indexinintvec++;
    std::vector<int> disjregs(numdisjregs); std::vector<int> numelems(numdisjregs);
    
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
        std::cout << "Error in 'rawfield' object: the mesh used to write file '" << filename << "' is not the same as the current one" << std::endl;
        abort();
    }
    
    
    ///// ALL CHECKS PERFORMED. LOAD THE DATA.
    
    // Calculate the number of sons:
    int numsons = (intdata.size() - indexinintvec) / (2*numdisjregs);
    // Get all the sons in this rawfield:
    std::vector<std::pair<std::vector<int>, shared_ptr<rawfield>>> allsons = getallsons();
    
    for (int i = 0; i < numsons; i++)
    {
        shared_ptr<rawfield> curson = allsons[i].second;
        shared_ptr<coefmanager> curcoefmanager = curson->mycoefmanager;
    
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
}
        

std::vector<std::vector<densematrix>> rawfield::interpolate(int whichderivative, int formfunctioncomponent, elementselector& elemselect, std::vector<double>& evaluationcoordinates)
{
    // Get all disjoint regions in the element selector:
    std::vector<int> alldisjregs = elemselect.getdisjointregions();
    
    std::vector<std::vector<densematrix>> out = {};

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
        
        int elementtypenumber = universe::mymesh->getdisjointregions()->getelementtypenumber(mydisjregs[0]);
        std::vector<std::vector<densematrix>> currentinterp = interpolate(whichderivative, formfunctioncomponent, elementtypenumber, elemselect.gettotalorientation(), getinterpolationorder(mydisjregs[0]), elemselect.getelementnumbers(), evaluationcoordinates);
        
        // Preallocate the harmonics not in 'out':
        if (out.size() < currentinterp.size())
            out.resize(currentinterp.size());
        for (int h = 0; h < currentinterp.size(); h++)
        {
            if (currentinterp[h].size() == 1 && out[h].size() == 0)
                out[h] = {densematrix(elemselect.countincurrentorientation(), evaluationcoordinates.size()/3, 0)};
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

densematrix rawfield::getcoefficients(int elementtypenumber, int interpolorder, std::vector<int> elementlist)
{   
    disjointregions* mydisjointregions = universe::mymesh->getdisjointregions();  
    elements* myelements = universe::mymesh->getelements();
    
    // In case the field is the x, y or z coordinate field:
    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {
        int mycoordinate = 0;
        if (mytypename == "y")
            mycoordinate = 1;
        if (mytypename == "z")
            mycoordinate = 2;
        
        std::vector<double>* mynodecoordinates = universe::mymesh->getnodes()->getcoordinates();
     
        element myelement(elementtypenumber, myelements->getcurvatureorder());        
        int numcurvednodes = myelement.countcurvednodes();
        
        int numrows = numcurvednodes;
        int numcols = elementlist.size();
        
        densematrix coefmatrix(numrows, numcols);
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

        densematrix coefmatrix(numrows, numcols);
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

std::vector<std::vector<densematrix>> rawfield::interpolate(int whichderivative, int formfunctioncomponent, int elementtypenumber, int totalorientation, int interpolorder, std::vector<int> elementnumbers, std::vector<double>& evaluationcoordinates)
{   
    elements* myelements = universe::mymesh->getelements();

    if (mytypename == "x" || mytypename == "y" || mytypename == "z")
    {        
        // Get the coefficients. The interpolation order is not used here.
        densematrix mycoefs = getcoefficients(elementtypenumber, -1, elementnumbers);
        // Compute the form functions evaluated at the evaluation points:
        lagrangeformfunction mylagrange(elementtypenumber, myelements->getcurvatureorder(), evaluationcoordinates);
        densematrix myformfunctionvalue = mylagrange.getderivative(whichderivative);

        mycoefs.transpose();
        densematrix coefstimesformfunctions = mycoefs.multiply(myformfunctionvalue);

        return {{},{coefstimesformfunctions}};
    }
    else
    {
        if (myharmonics.size() == 0)
        {
            // Get the coefficients:
            densematrix mycoefs = getcoefficients(elementtypenumber, interpolorder, elementnumbers);
            // Compute the form functions evaluated at the evaluation points.
            // This reuses as much as possible what's already been computed:
            hierarchicalformfunctioncontainer* val = universe::gethff(mytypename, elementtypenumber, interpolorder, evaluationcoordinates);
            
            // We can get the total orientation from elementnumbers[0] since
            // we require that all elements have the same total orientation.
            densematrix myformfunctionvalue = val->tomatrix(totalorientation, interpolorder, whichderivative, formfunctioncomponent);
            
            mycoefs.transpose();
            densematrix coefstimesformfunctions = mycoefs.multiply(myformfunctionvalue);

            return {{},{coefstimesformfunctions}};
        }
        else
        {
            std::vector<std::vector<densematrix>> coefstimesformfunctions(myharmonics.size());

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


