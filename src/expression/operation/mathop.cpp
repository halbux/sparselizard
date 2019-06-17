#include "mathop.h"


double mathop::getpi(void)
{
    return 3.1415926535897932384;
}

int mathop::regionunion(const std::vector<int> physregs)
{
    return (universe::mymesh->getphysicalregions())->createunion(physregs);
}

int mathop::regionintersection(const std::vector<int> physregs)
{
    return (universe::mymesh->getphysicalregions())->createintersection(physregs);
}

int mathop::regionexclusion(int physreg, int toexclude)
{
    return (universe::mymesh->getphysicalregions())->createexclusion(physreg, toexclude);
}

void mathop::printvector(std::vector<double> input)
{
    std::cout << "Vector size is " << input.size() << std::endl;
    for (int i = 0; i < input.size(); i++)
        std::cout << input[i] << " ";
    std::cout << std::endl;
}

void mathop::printvector(std::vector<int> input)
{
    std::cout << "Vector size is " << input.size() << std::endl;
    for (int i = 0; i < input.size(); i++)
        std::cout << input[i] << " ";
    std::cout << std::endl;
}

void mathop::writevector(std::string filename, std::vector<double>& towrite)
{
	if (towrite.size() == 0)
		return;

	// 'file' cannot take a std::string argument --> filename.c_str():
	std::ofstream name (filename.c_str());
	if (name.is_open())
	{
		// To write all doubles with enough digits to the file:
		name << std::setprecision(17);
		
		for (int i = 0; i < towrite.size()-1; i++)
			name << towrite[i] << ",";
		name << towrite[towrite.size()-1];
		
		name.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << filename << " or file not found" << std::endl;
		abort();
	}
}

expression mathop::norm(expression expr)
{
    if (expr.countcolumns() > 1)
    {
        std::cout << "Error in 'mathop' namespace: can only compute the norm of column vectors" << std::endl;
        abort();
    }

    expression mynorm = pow(expr.at(0,0),2);
    for (int i = 1; i < expr.countrows(); i++)
		mynorm = mynorm + pow(expr.at(i,0),2);

	return sqrt(mynorm);
}

expression mathop::normal(int surfphysreg)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    int elementdimension = universe::mymesh->getphysicalregions()->get(surfphysreg)->getelementdimension();
    
    if (problemdimension-1 != elementdimension || problemdimension == 1)
    {
        std::cout << "Error in 'mathop' namespace: can only compute the normal to a surface in 3D and to a line in 2D" << std::endl;
        abort();
    }
    
    expression expr;
    if (problemdimension == 2)
    {
        expression mynorm = sqrt(expr.invjac(0,1)*expr.invjac(0,1)+expr.invjac(1,1)*expr.invjac(1,1));
        mynorm.reuseit();
		if (universe::isaxisymmetric)
        	return array3x1(expr.invjac(0,1), expr.invjac(1,1), 0)/mynorm;
		else
			return array2x1(expr.invjac(0,1), expr.invjac(1,1))/mynorm;
    }
    if (problemdimension == 3)
    {
        expression mynorm = sqrt(expr.invjac(0,2)*expr.invjac(0,2)+expr.invjac(1,2)*expr.invjac(1,2)+expr.invjac(2,2)*expr.invjac(2,2));
        mynorm.reuseit();
        return array3x1(expr.invjac(0,2), expr.invjac(1,2), expr.invjac(2,2))/mynorm;
    }
}

void mathop::scatterwrite(std::string filename, std::vector<double> xcoords, std::vector<double> ycoords, std::vector<double> zcoords, std::vector<double> compxevals, std::vector<double> compyevals, std::vector<double> compzevals)
{
    int n = xcoords.size();

    if (n == 0)
        return;

    // Is the data to write scalar or a vector:
    bool isscalar = (compxevals.size() > 0 && compyevals.size() == 0 && compzevals.size() == 0);

    if (isscalar == false && compxevals.size() == 0)
        compxevals = std::vector<double>(n,0);
    if (isscalar == false && compyevals.size() == 0)
        compyevals = std::vector<double>(n,0);
    if (isscalar == false && compzevals.size() == 0)
        compzevals = std::vector<double>(n,0);

    if (xcoords.size() != n || ycoords.size() != n || zcoords.size() != n || compxevals.size() != n || (isscalar == false && (compyevals.size() != n || compzevals.size() != n)))
    {
        std::cout << "Error in 'mathop' namespace: size of 'scatterwrite' arguments do not match" << std::endl;
        abort();
    }
    
    iodata datatowrite(1, 1, isscalar, {});
    datatowrite.addcoordinates(0, densematrix(n,1,xcoords), densematrix(n,1,ycoords), densematrix(n,1,zcoords));
    if (isscalar)
    	datatowrite.adddata(0, {densematrix(n,1,compxevals)});
	else
		datatowrite.adddata(0, {densematrix(n,1,compxevals), densematrix(n,1,compyevals), densematrix(n,1,compzevals)});
		
	iointerface::writetofile(filename, datatowrite);
}
    
    
void mathop::setaxisymmetry(void) { universe::isaxisymmetric = true; }

void mathop::setfundamentalfrequency(double f) { universe::fundamentalfrequency = f; }
void mathop::settime(double t) { universe::currenttimestep = t; }
double mathop::gettime(void) { return universe::currenttimestep; }

expression mathop::t(void) { expression exp; return exp.time(); }

field mathop::elementsize(int physreg)
{
	int problemdimension = universe::mymesh->getmeshdimension();
	int elementdimension = universe::mymesh->getphysicalregions()->get(physreg)->getelementdimension();

	if (problemdimension != elementdimension)
	{
		std::cout << "Error in 'mathop' namespace: trying to get the size of a " << elementdimension << "D element in a " << problemdimension << "D problem (dimensions must be equal)" << std::endl;
		abort();
	}

	// Create a 'one' type field that stores the element volume/area/length:
	field one("one");
	
	// Create a formulation to integrate a constant one and get the size:
	formulation formul;
	formul += integration(physreg, -tf(one));
	formul.generate();
	
	vec rhs = formul.rhs();
	
	// Save to the one field:
	one.setdata(physreg, rhs);
	
	return one;
}

expression mathop::dx(expression input) { return input.spacederivative(1); }
expression mathop::dy(expression input) { return input.spacederivative(2); }
expression mathop::dz(expression input) { return input.spacederivative(3); }

expression mathop::dt(expression input) { return input.timederivative(1); }
expression mathop::dtdt(expression input) { return input.timederivative(2); }

expression mathop::sin(expression input) { return input.sin(); }
expression mathop::cos(expression input) { return input.cos(); }
expression mathop::tan(expression input) { return input.tan(); }
expression mathop::asin(expression input) { return input.asin(); }
expression mathop::acos(expression input) { return input.acos(); }
expression mathop::atan(expression input) { return input.atan(); }
expression mathop::abs(expression input) { return input.abs(); }
expression mathop::sqrt(expression input) { return pow(input, 0.5); }
expression mathop::log10(expression input) { return input.log10(); }
expression mathop::pow(expression base, expression exponent) { return base.pow(exponent); }

expression mathop::comp(int selectedcomp, expression input) 
{ 
    std::vector<expression> mycomp(input.countcolumns());
    for (int i = 0; i < input.countcolumns(); i++)
        mycomp[i] = input.at(selectedcomp,i);
    return expression(1, input.countcolumns(), mycomp);
}

expression mathop::compx(expression input) { return comp(0,input); }
expression mathop::compy(expression input) { return comp(1,input); }
expression mathop::compz(expression input) { return comp(2,input); }

expression mathop::entry(int row, int col, expression input) { return input.at(row,col); }

expression mathop::transpose(expression input) { return input.transpose(); }
expression mathop::inverse(expression input) { return input.invert(); }
expression mathop::determinant(expression input) { return input.determinant(); }

expression mathop::grad(expression input)
{
    if (input.countcolumns() != 1 || input.countrows() > 3)
    {
        std::cout << "Error in 'mathop' namespace: can only take the gradient of a scalar or an up to length 3 column vector" << std::endl;
        abort();
    }
    
    // Cylindrical transformation of the gradient of a vector (different than of a scalar):
    if (universe::isaxisymmetric && input.countrows() > 1)
    {
    	field x("x");
    	// The output must be a nx3 matrix:
    	if (input.countrows() == 2)
    		return expression(2,3, {compx(dx(input)),compx(dy(input)),0, compy(dx(input)),compy(dy(input)),0});
		if (input.countrows() == 3)
			return expression(3,3, {compx(dx(input)),compx(dy(input)),-1.0/x*compz(input), compy(dx(input)),compy(dy(input)),0, compz(dx(input)),compz(dy(input)),1.0/x*compx(input)});
    }

    int problemdimension = universe::mymesh->getmeshdimension();
    // In case of axisymmetry we need a 3 component output vector:
	if (universe::isaxisymmetric)
		problemdimension++;

    std::vector<expression> myexprs = {};
    for (int comp = 0; comp < input.countrows(); comp++)
    {
        for (int i = 0; i < problemdimension; i++)
            myexprs.push_back(input.spacederivative(i+1).at(comp,0));
    }
    
    expression output(input.countrows(), problemdimension, myexprs);
    
    // We want the gradient of a scalar to be a column vector:
    if (input.countrows() == 1)
    	output = output.transpose();
    
    return output;
}

expression mathop::div(expression input)
{
    if (input.countcolumns() != 1 || input.countrows() > 3)
    {
        std::cout << "Error in 'mathop' namespace: can only take the divergence of an up to length 3 column vector" << std::endl;
        abort();
    }
    
    if (universe::isaxisymmetric)
    {
    	field x("x");
    	switch (input.countrows())
    	{
			case 1:
				return 1.0/x*compx(input)+compx(dx(input));
			case 2:
				return 1.0/x*compx(input)+compx(dx(input))+compy(dy(input));
			case 3:
				return 1.0/x*compx(input)+compx(dx(input))+compy(dy(input));
    	}
	}

	switch (input.countrows())
	{
		case 1:
		    return dx(input);
		case 2:
		    return compx(dx(input))+compy(dy(input));
		case 3:
		    return compx(dx(input))+compy(dy(input))+compz(dz(input));
	}

}

expression mathop::curl(expression input)
{
    bool ishcurlfield = input.isprojectedfield();
    input = input.getunprojectedfield();
    
    if (input.countcolumns() > 1 || input.countrows() > 3)
    {
        std::cout << "Error in 'mathop' namespace: can only take the curl of an up to length 3 column vector" << std::endl;
        abort();
    }

    // The curl of a hcurl type field is computed in a special way:
    if (ishcurlfield == false)
    {
    	if (universe::isaxisymmetric)
    	{
    		field x("x");
		    switch (input.countrows())
		    {
		        case 1:
		            return expression(3,1,{0, 0, compx(dy(input))});
		        case 2:
		            return expression(3,1,{0, 0, compx(dy(input))-compy(dx(input))});
		        case 3:
		            return expression(3,1,{compz(dy(input)), -1.0/x*compz(input)-compz(dx(input)), -compx(dy(input))+compy(dx(input))});
		    }
    	}
    	
        switch (input.countrows())
        {
            case 1:
                return expression(3,1,{0, 0, 0});
            case 2:
                return expression(3,1,{0, 0, compy(dx(input))-compx(dy(input))});
            case 3:
                return expression(3,1,{compz(dy(input))-compy(dz(input)), compx(dz(input))-compz(dx(input)), compy(dx(input))-compx(dy(input))});
        }
    }
    else
    {
        expression expr;
        // We should always have 3 components.
        // This is the curl in the reference element:
        expr = expression(3,1,{compz(input.kietaphiderivative(2))-compy(input.kietaphiderivative(3)), compx(input.kietaphiderivative(3))-compz(input.kietaphiderivative(1)), compy(input.kietaphiderivative(1))-compx(input.kietaphiderivative(2))});

        // The curl of a 1 form (i.e. hcurl type field) is brought back 
        // from the reference element with the following transformation:
        expr = transpose(expr.jac())*expr/expr.detjac();
        return expr;
    }
}

expression mathop::crossproduct(expression a, expression b)
{
    if (a.countcolumns() != 1 || b.countcolumns() != 1 || a.countrows() > 3 || b.countrows() > 3)
    {
        std::cout << "Error in 'mathop' namespace: can only take the cross product of up to length 3 column vectors" << std::endl;
        abort();
    }

	a = a.resize(3,1);
	b = b.resize(3,1);

	expression a1 = compx(a), a2 = compy(a), a3 = compz(a);
	expression b1 = compx(b), b2 = compy(b), b3 = compz(b);

	expression crossprodexpr = expression(3,1, { a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1 });

	return crossprodexpr;
}

expression mathop::frobeniusproduct(expression a, expression b)
{
    if (a.countcolumns() != b.countcolumns() || a.countrows() != b.countrows())
    {
        std::cout << "Error in 'mathop' namespace: dimension mismatch for Frobenius product" << std::endl;
        abort();
    }
    
    expression output;
    for (int i = 0; i < a.countrows(); i++)
    {
        for (int j = 0; j < a.countcolumns(); j++)
        {
            if (i == 0 && j == 0)
                output = a.at(i,j) * b.at(i,j);
            else
                output = output + a.at(i,j) * b.at(i,j);
        }
    }

    return output;
}

expression mathop::trace(expression a)
{
    if (a.countcolumns() != a.countrows())
    {
        std::cout << "Error in 'mathop' namespace: can only get the trace of a square matrix" << std::endl;
        abort();
    }
    
    expression output;
    for (int i = 0; i < a.countrows(); i++)
    {
        if (i == 0)
            output = a.at(i,i);
        else
            output = output + a.at(i,i);
    }

    return output;
}

expression mathop::detjac(void) 
{ 
	expression expr;
	return expr.detjac(); 
}

integration mathop::integral(int physreg, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    return integration(physreg, tointegrate, integrationorderdelta, blocknumber);
}

integration mathop::integral(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    return integration(physreg, meshdeform, tointegrate, integrationorderdelta, blocknumber);
}

integration mathop::integral(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    return integration(physreg, numcoefharms, tointegrate, integrationorderdelta, blocknumber);
}

integration mathop::integral(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    return integration(physreg, numcoefharms, meshdeform, tointegrate, integrationorderdelta, blocknumber);
}

expression mathop::dof(expression input, int physreg) { return input.dof(physreg); }
expression mathop::tf(expression input, int physreg) { return input.tf(physreg); }



expression mathop::array1x1(expression term11) 
{
	std::vector<expression> terms = {term11};
    return expression(1,1, terms);
}

expression mathop::array1x2(expression term11, expression term12) 
{
    return expression(1,2, {term11, term12});
}

expression mathop::array1x3(expression term11, expression term12, expression term13) 
{
    return expression(1,3, {term11, term12, term13});
}

expression mathop::array2x1(expression term11, expression term21) 
{
    return expression(2,1, {term11, term21});
}

expression mathop::array2x2(expression term11, expression term12, expression term21, expression term22) 
{
    return expression(2,2, {term11, term12, term21, term22});
}

expression mathop::array2x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23) 
{
    return expression(2,3, {term11, term12, term13, term21, term22, term23});
}

expression mathop::array3x1(expression term11, expression term21, expression term31) 
{
    return expression(3,1, {term11, term21, term31});
}

expression mathop::array3x2(expression term11, expression term12, expression term21, expression term22, expression term31, expression term32) 
{
    return expression(3,2, {term11, term12, term21, term22, term31, term32});
}

expression mathop::array3x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23, expression term31, expression term32, expression term33) 
{
    return expression(3,3, {term11, term12, term13, term21, term22, term23, term31, term32, term33});
}
 
vec mathop::solve(mat A, vec b, std::string soltype, bool diagscaling)
{
    if (soltype != "lu")
    {
        std::cout << "Error in 'mathop' namespace: unknown direct solver type '" << soltype << "' (use 'lu')" << std::endl;
        abort();
    }
    if (A.countrows() != b.size())
    {
        std::cout << "Error in 'mathop' namespace: direct solve of Ax = b failed (size of A and b do not match)" << std::endl;
        abort();
    }
    
    if (A.getpointer() == NULL || b.getpointer() == NULL)
    {
        std::cout << "Error in 'mathop' namespace: direct solve of Ax = b failed (A or b is undefined)" << std::endl;
        abort();
    }
    
    // The copy of the rhs is returned in case there is no nonzero entry in A:
    if (A.countnnz() == 0)
    	return b.copy();

    Vec bpetsc = b.getpetsc();
    Mat Apetsc = A.getpetsc();

    vec sol(shared_ptr<rawvec>(new rawvec(b.getpointer()->getdofmanager())));
    Vec solpetsc = sol.getpetsc();

    KSP* ksp = A.getpointer()->getksp();
    
    if (A.getpointer()->isludefined() == false)
    {
        PC pc;
        KSPCreate(PETSC_COMM_WORLD, ksp);
        KSPSetOperators(*ksp, Apetsc, Apetsc);
        // Perform a diagonal scaling for improved matrix conditionning.
        // This modifies the matrix A and right handside b!
        if (diagscaling == true)
            KSPSetDiagonalScale(*ksp, PETSC_TRUE);
        KSPSetFromOptions(*ksp);

        KSPGetPC(*ksp,&pc);
        PCSetType(pc,PCLU);
        PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
    }
    	
    KSPSolve(*ksp, bpetsc, solpetsc);
   
    A.getpointer()->isludefined(true);
    
    if (A.getpointer()->islutobereused() == false)
    {
        KSPDestroy(ksp);
        A.getpointer()->isludefined(false);
    }
    
    return sol;
}

int mykspmonitor(KSP ksp, PetscInt iter, PetscReal resnorm, void* unused)
{
    std::cout << iter << " KSP residual norm " << resnorm << std::endl;
    return 0;
}

void mathop::solve(mat A, vec b, vec sol, double& relrestol, int& maxnumit, std::string soltype, std::string precondtype, int verbosity, bool diagscaling)
{
    if (soltype != "gmres" && soltype != "bicgstab")
    {
        std::cout << "Error in 'mathop' namespace: unknown iterative solver type '" << soltype << "' (use 'gmres' or 'bicgstab')" << std::endl;
        abort();
    }
    if (precondtype != "ilu" && precondtype != "sor")
    {
        std::cout << "Error in 'mathop' namespace: unknown preconditioner type '" << precondtype << "' (use 'ilu' or 'sor')" << std::endl;
        abort();
    }
    if (A.countrows() != b.size())
    {
        std::cout << "Error in 'mathop' namespace: iterative solve of Ax = b failed (size of A and b do not match)" << std::endl;
        abort();
    }
    if (A.countrows() != sol.size())
    {
        std::cout << "Error in 'mathop' namespace: iterative solve of Ax = b failed (size of A and x do not match)" << std::endl;
        abort();
    }
    
    if (A.getpointer() == NULL || b.getpointer() == NULL || sol.getpointer() == NULL)
    {
        std::cout << "Error in 'mathop' namespace: iterative solve of Ax = b failed (A, x or b is undefined)" << std::endl;
        abort();
    }

    Vec bpetsc = b.getpetsc();
    Mat Apetsc = A.getpetsc();

    Vec solpetsc = sol.getpetsc();

    KSP* ksp = A.getpointer()->getksp();
    
    KSPCreate(PETSC_COMM_WORLD, ksp);
    KSPSetOperators(*ksp, Apetsc, Apetsc);
    // Perform a diagonal scaling for improved matrix conditionning.
    // This modifies the matrix A and right handside b!
    if (diagscaling == true)
        KSPSetDiagonalScale(*ksp, PETSC_TRUE);
    
    if (soltype == "gmres")
        KSPSetType(*ksp, KSPGMRES);
    if (soltype == "bicgstab")
        KSPSetType(*ksp, KSPBCGS);
    
    // The initial guess is provided in vector sol:
    KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE);
        
    KSPSetTolerances(*ksp, relrestol, PETSC_DEFAULT, PETSC_DEFAULT, maxnumit);
        
    // Request to print the iteration information to the console:
    if (verbosity > 0)
        KSPMonitorSet(*ksp, mykspmonitor, PETSC_NULL, PETSC_NULL);
        
    KSPSetFromOptions(*ksp);
    	
    // Use a preconditioner:
    PC pc;
    KSPGetPC(*ksp,&pc);
    if (precondtype == "ilu")
        PCSetType(pc,PCILU);
    if (precondtype == "sor")
        PCSetType(pc,PCSOR);
    	
    KSPSolve(*ksp, bpetsc, solpetsc);
    
    // Get the number of required iterations and the residual norm:
    KSPGetIterationNumber(*ksp, &maxnumit);
    KSPGetResidualNorm(*ksp, &relrestol);
   
    KSPDestroy(ksp);
}

void mathop::solve(formulation formul)
{
	// Get all fields in the formulation:
	std::vector<shared_ptr<rawfield>> allfields = formul.getdofmanager()->getfields();

	// Remove leftovers (if any):
	mat A = formul.A(); vec b = formul.b();
	// Generate:
	formul.generate();
	// Solve:
	vec sol = mathop::solve(formul.A(), formul.b());
	
	// Save to fields:
	for (int i = 0; i < allfields.size(); i++)
		allfields[i]->setdata(-1, sol|field(allfields[i]));
}

void mathop::solve(std::vector<formulation> formuls)
{
	for (int i = 0; i < formuls.size(); i++)
		solve(formuls[i]);
}



////////// PREDEFINED OPERATORS

expression mathop::strain(expression input)
{
    if ((input.countrows() != 2 && input.countrows() != 3) || input.countcolumns() != 1)
    {
        std::cout << "Error in 'mathop' namespace: can only compute the strains of a 2x1 or 3x1 column vector" << std::endl;
        abort();
    }
    
    expression gradu = mathop::grad(input);
    
	if (input.countrows() == 2)
		return expression(3,1,{gradu.at(0,0), gradu.at(1,1), gradu.at(1,0) + gradu.at(0,1)});
	if (input.countrows() == 3)
		return expression(6,1,{gradu.at(0,0), gradu.at(1,1), gradu.at(2,2), gradu.at(2,1) + gradu.at(1,2), gradu.at(0,2) + gradu.at(2,0), gradu.at(0,1) + gradu.at(1,0)});
}

expression mathop::greenlagrangestrain(expression gradu)
{
	// This can be called since gradu is nonlinear in u and can thus not include a dof or tf:
	gradu.reuseit();

	if (gradu.countrows() == 2 && gradu.countcolumns() == 2)
	{
		expression dxcompxu = entry(0,0,gradu), dxcompyu = entry(1,0,gradu);
		expression dycompxu = entry(0,1,gradu), dycompyu = entry(1,1,gradu);

		expression output = expression(3,1, {
								dxcompxu + 0.5*(pow(dxcompxu,2) + pow(dxcompyu,2)), 
								dycompyu + 0.5*(pow(dycompxu,2) + pow(dycompyu,2)), 
								dycompxu + dxcompyu + dxcompxu * dycompxu + dxcompyu * dycompyu});
		output.reuseit();
		return output;
	}
	if (gradu.countrows() == 3 && gradu.countcolumns() == 3)
	{
		expression dxcompxu = entry(0,0,gradu), dxcompyu = entry(1,0,gradu), dxcompzu = entry(2,0,gradu);
		expression dycompxu = entry(0,1,gradu), dycompyu = entry(1,1,gradu), dycompzu = entry(2,1,gradu);
		expression dzcompxu = entry(0,2,gradu), dzcompyu = entry(1,2,gradu), dzcompzu = entry(2,2,gradu);

		expression output = expression(6,1, {
								dxcompxu + 0.5*(pow(dxcompxu,2) + pow(dxcompyu,2) + pow(dxcompzu,2)), 
								dycompyu + 0.5*(pow(dycompxu,2) + pow(dycompyu,2) + pow(dycompzu,2)), 
								dzcompzu + 0.5*(pow(dzcompxu,2) + pow(dzcompyu,2) + pow(dzcompzu,2)), 
								dzcompyu + dycompzu + dycompxu * dzcompxu + dycompyu * dzcompyu + dycompzu * dzcompzu, 
								dzcompxu + dxcompzu + dxcompxu * dzcompxu + dxcompyu * dzcompyu + dxcompzu * dzcompzu, 
								dycompxu + dxcompyu + dxcompxu * dycompxu + dxcompyu * dycompyu + dxcompzu * dycompzu});
		output.reuseit();
		return output;
	}

    std::cout << "Error in 'mathop' namespace: expected a 2x2 or 3x3 matrix as input for greenlagrangestrain()" << std::endl;
    abort();
}

expression mathop::vonmises(expression stress)
{
	if (stress.countcolumns() != 1 || stress.countrows() != 6)
	{
		std::cout << "Error in 'mathop' namespace: expected the 3D stress tensor in Voigt notation (6 rows, 1 column)" << std::endl;
		abort();
	}

	expression s11 = stress.at(0,0), s22 = stress.at(1,0), s33 = stress.at(2,0), s23 = stress.at(3,0), s13 = stress.at(4,0), s12 = stress.at(5,0);
	s11.reuseit(); s22.reuseit(); s33.reuseit();

	return sqrt( 0.5*( pow(s11-s22,2)+pow(s22-s33,2)+pow(s33-s11,2) ) + 3.0*( pow(s12,2)+pow(s23,2)+pow(s13,2) ) );
}

expression mathop::predefinedmassconservation(expression dofv, expression tfp, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant)
{
    if (isdensityconstant)
        return div(dofv)*tfp;
        
    if (includetimederivs)
        return ( rho*div(dofv)*tfp + dofv*gradrho*tfp + dtrho*tfp );
    else
        return ( rho*div(dofv)*tfp + dofv*gradrho*tfp );
}

expression mathop::predefinedinertialforce(expression dofv, expression tfv, expression v, expression rho)
{
    rho.reuseit(); v.reuseit();
    
    return ( -rho*( grad(v)*dofv + grad(dofv)*v - grad(v)*v )*tfv );
}

expression mathop::predefinedviscousforce(expression dofv, expression tfv, expression mu, bool isdensityconstant, bool isviscosityconstant)
{
    mu.reuseit();
    
    if (isdensityconstant && isviscosityconstant)
        return ( - mu*frobeniusproduct(grad(dofv), grad(tfv)) );
    
    if (isdensityconstant)
        return ( - mu*frobeniusproduct(grad(dofv), grad(tfv)) - mu*frobeniusproduct(transpose(grad(dofv)), grad(tfv)) );
    else
        return ( - mu*frobeniusproduct(grad(dofv), grad(tfv)) - mu*frobeniusproduct(transpose(grad(dofv)), grad(tfv)) + (2.0/3.0)*mu*div(dofv)*trace(grad(tfv)) );
}


////////// PREDEFINED FORMULATIONS

expression mathop::predefinedelasticity(expression dofu, expression tfu, expression E, expression nu, std::string myoption)
{
	// Hooke's matrix:
	expression H(6,6, {1-nu,nu,nu,0,0,0,  nu,1-nu,nu,0,0,0,  nu,nu,1-nu,0,0,0,  0,0,0,0.5*(1-2*nu),0,0,  0,0,0,0,0.5*(1-2*nu),0,  0,0,0,0,0,0.5*(1-2*nu)});
	expression coef = E/(1+nu)/(1-2*nu);
	coef.reuseit();
	H = coef * H;
	return predefinedelasticity(dofu, tfu, H, myoption);
}

expression mathop::predefinedelasticity(expression dofu, expression tfu, expression H, std::string myoption)
{
    if (dofu.countrows() != tfu.countrows() || dofu.countcolumns() != 1 || tfu.countcolumns() != 1 || dofu.countrows() == 1)
    {
        std::cout << "Error in 'mathop' namespace: first arguments in 'predefinedelasticity' must be either 2x1 or 3x1 vectors" << std::endl;
        abort();
    }
    if (dofu.countrows() == 2)
	{
		if (myoption == "planestrain")
			H = expression(3,3, {H.at(0,0),H.at(0,1),H.at(0,5), H.at(1,0),H.at(1,1),H.at(1,5), H.at(5,0),H.at(5,1),H.at(5,5) });
		if (myoption == "planestress")
		{
			expression subdet,ezztoexx,ezztoeyy,ezztog12,g23toexx,g23toeyy,g23tog12,g13toexx,g13toeyy,g13tog12;

			subdet = H.at(2,2)*H.at(3,3)*H.at(4,4)-H.at(2,2)*H.at(3,4)*H.at(4,3)-H.at(2,3)*H.at(3,2)*H.at(4,4)+H.at(2,3)*H.at(3,4)*H.at(4,2)+H.at(2,4)*H.at(3,2)*H.at(4,3)-H.at(2,4)*H.at(3,3)*H.at(4,2);
			subdet.reuseit();

			// This is the extra contribution of ezz: 
			ezztoexx = ( H.at(3,0)*( H.at(2,3)*H.at(4,4) - H.at(2,4)*H.at(4,3)) - H.at(4,0)*( H.at(2,3)*H.at(3,4) - H.at(2,4)*H.at(3,3)) - H.at(2,0)*( H.at(3,3)*H.at(4,4) - H.at(3,4)*H.at(4,3)))/subdet;
			ezztoeyy = ( H.at(3,1)*( H.at(2,3)*H.at(4,4) - H.at(2,4)*H.at(4,3)) - H.at(4,1)*( H.at(2,3)*H.at(3,4) - H.at(2,4)*H.at(3,3)) - H.at(2,1)*( H.at(3,3)*H.at(4,4) - H.at(3,4)*H.at(4,3)))/subdet;
			ezztog12 = ( H.at(3,5)*( H.at(2,3)*H.at(4,4) - H.at(2,4)*H.at(4,3)) - H.at(4,5)*( H.at(2,3)*H.at(3,4) - H.at(2,4)*H.at(3,3)) - H.at(2,5)*( H.at(3,3)*H.at(4,4) - H.at(3,4)*H.at(4,3)))/subdet;
			// This is the extra contribution of g23:
			g23toexx = ( H.at(4,0)*( H.at(2,2)*H.at(3,4) - H.at(2,4)*H.at(3,2)) - H.at(3,0)*( H.at(2,2)*H.at(4,4) - H.at(2,4)*H.at(4,2)) + H.at(2,0)*( H.at(3,2)*H.at(4,4) - H.at(3,4)*H.at(4,2)))/subdet;
			g23toeyy = ( H.at(4,1)*( H.at(2,2)*H.at(3,4) - H.at(2,4)*H.at(3,2)) - H.at(3,1)*( H.at(2,2)*H.at(4,4) - H.at(2,4)*H.at(4,2)) + H.at(2,1)*( H.at(3,2)*H.at(4,4) - H.at(3,4)*H.at(4,2)))/subdet;
			g23tog12 = ( H.at(4,5)*( H.at(2,2)*H.at(3,4) - H.at(2,4)*H.at(3,2)) - H.at(3,5)*( H.at(2,2)*H.at(4,4) - H.at(2,4)*H.at(4,2)) + H.at(2,5)*( H.at(3,2)*H.at(4,4) - H.at(3,4)*H.at(4,2)))/subdet;
			// This is the extra contribution of g13: 
			g13toexx = ( H.at(3,0)*( H.at(2,2)*H.at(4,3) - H.at(2,3)*H.at(4,2)) - H.at(4,0)*( H.at(2,2)*H.at(3,3) - H.at(2,3)*H.at(3,2)) - H.at(2,0)*( H.at(3,2)*H.at(4,3) - H.at(3,3)*H.at(4,2)))/subdet;
			g13toeyy = ( H.at(3,1)*( H.at(2,2)*H.at(4,3) - H.at(2,3)*H.at(4,2)) - H.at(4,1)*( H.at(2,2)*H.at(3,3) - H.at(2,3)*H.at(3,2)) - H.at(2,1)*( H.at(3,2)*H.at(4,3) - H.at(3,3)*H.at(4,2)))/subdet;
			g13tog12 = ( H.at(3,5)*( H.at(2,2)*H.at(4,3) - H.at(2,3)*H.at(4,2)) - H.at(4,5)*( H.at(2,2)*H.at(3,3) - H.at(2,3)*H.at(3,2)) - H.at(2,5)*( H.at(3,2)*H.at(4,3) - H.at(3,3)*H.at(4,2)))/subdet;

			ezztoexx.reuseit(); ezztoeyy.reuseit(); ezztog12.reuseit(); g23toexx.reuseit(); g23toeyy.reuseit(); g23tog12.reuseit(); g13toexx.reuseit(); g13toeyy.reuseit(); g13tog12.reuseit();

			H = expression(3,3,{  
			H.at(0,0) + H.at(0,2)*ezztoexx+H.at(0,3)*g23toexx+H.at(0,4)*g13toexx,
			H.at(0,1) + H.at(0,2)*ezztoeyy+H.at(0,3)*g23toeyy+H.at(0,4)*g13toeyy,
			H.at(0,5) + H.at(0,2)*ezztog12+H.at(0,3)*g23tog12+H.at(0,4)*g13tog12,

			H.at(1,0) + H.at(1,2)*ezztoexx+H.at(1,3)*g23toexx+H.at(1,4)*g13toexx,
			H.at(1,1) + H.at(1,2)*ezztoeyy+H.at(1,3)*g23toeyy+H.at(1,4)*g13toeyy,
			H.at(1,5) + H.at(1,2)*ezztog12+H.at(1,3)*g23tog12+H.at(1,4)*g13tog12,

			H.at(5,0) + H.at(5,2)*ezztoexx+H.at(5,3)*g23toexx+H.at(5,4)*g13toexx,
			H.at(5,1) + H.at(5,2)*ezztoeyy+H.at(5,3)*g23toeyy+H.at(5,4)*g13toeyy,
			H.at(5,5) + H.at(5,2)*ezztog12+H.at(5,3)*g23tog12+H.at(5,4)*g13tog12
			});
		}
		
		if (myoption == "planestrain" || myoption == "planestress")
			return -( H *strain(dofu) )*strain(tfu);

		// If the option is not valid:
		std::cout << "Error in 'mathop' namespace: invalid option or no option provided for the 2D problem in 'predefinedelasticity'" << std::endl;
		std::cout << "Available choices are: 'planestrain', 'planestress'" << std::endl;
		abort();
	}
    if (dofu.countrows() == 3)
    {
        if (myoption.length() > 0)
        {
            std::cout << "Error in 'mathop' namespace: for a 3D problem the last string argument must be empty in 'predefinedelasticity'" << std::endl;
            abort();
        }
        return -( H*strain(dofu) ) * strain(tfu);
    }
}

expression mathop::predefinedelasticity(expression dofu, expression tfu, field u, expression E, expression nu, expression prestress, std::string myoption)
{
	// Hooke's matrix:
	expression H(6,6, {1-nu,nu,nu,0,0,0,  nu,1-nu,nu,0,0,0,  nu,nu,1-nu,0,0,0,  0,0,0,0.5*(1-2*nu),0,0,  0,0,0,0,0.5*(1-2*nu),0,  0,0,0,0,0,0.5*(1-2*nu)});
	expression coef = E/(1+nu)/(1-2*nu);
	coef.reuseit();
	H = coef * H;
	return predefinedelasticity(dofu, tfu, u, H, prestress, myoption);
}

expression mathop::predefinedelasticity(expression dofu, expression tfu, field u, expression H, expression prestress, std::string myoption)
{
    if (dofu.countrows() != tfu.countrows() || dofu.countcolumns() != 1 || tfu.countcolumns() != 1 || dofu.countrows() == 1)
    {
        std::cout << "Error in 'mathop' namespace: first arguments in 'predefinedelasticity' must be either 2x1 or 3x1 vectors" << std::endl;
        abort();
    }
    if (dofu.countrows() == 2)
	{
		if (prestress.iszero() == false && (prestress.countcolumns() != 1 || prestress.countrows() != 3))
		{
		    std::cout << "Error in 'mathop' namespace: expected a 3x1 sized prestress vector (Voigt form) in 'predefinedelasticity' (set scalar 0.0 if no prestress)" << std::endl;
		    abort();
		}

		if (myoption == "planestrain")
			H = expression(3,3, {H.at(0,0),H.at(0,1),H.at(0,5), H.at(1,0),H.at(1,1),H.at(1,5), H.at(5,0),H.at(5,1),H.at(5,5) });
		if (myoption == "planestress")
		{
			expression subdet,ezztoexx,ezztoeyy,ezztog12,g23toexx,g23toeyy,g23tog12,g13toexx,g13toeyy,g13tog12;

			subdet = H.at(2,2)*H.at(3,3)*H.at(4,4)-H.at(2,2)*H.at(3,4)*H.at(4,3)-H.at(2,3)*H.at(3,2)*H.at(4,4)+H.at(2,3)*H.at(3,4)*H.at(4,2)+H.at(2,4)*H.at(3,2)*H.at(4,3)-H.at(2,4)*H.at(3,3)*H.at(4,2);
			subdet.reuseit();

			// This is the extra contribution of ezz: 
			ezztoexx = ( H.at(3,0)*( H.at(2,3)*H.at(4,4) - H.at(2,4)*H.at(4,3)) - H.at(4,0)*( H.at(2,3)*H.at(3,4) - H.at(2,4)*H.at(3,3)) - H.at(2,0)*( H.at(3,3)*H.at(4,4) - H.at(3,4)*H.at(4,3)))/subdet;
			ezztoeyy = ( H.at(3,1)*( H.at(2,3)*H.at(4,4) - H.at(2,4)*H.at(4,3)) - H.at(4,1)*( H.at(2,3)*H.at(3,4) - H.at(2,4)*H.at(3,3)) - H.at(2,1)*( H.at(3,3)*H.at(4,4) - H.at(3,4)*H.at(4,3)))/subdet;
			ezztog12 = ( H.at(3,5)*( H.at(2,3)*H.at(4,4) - H.at(2,4)*H.at(4,3)) - H.at(4,5)*( H.at(2,3)*H.at(3,4) - H.at(2,4)*H.at(3,3)) - H.at(2,5)*( H.at(3,3)*H.at(4,4) - H.at(3,4)*H.at(4,3)))/subdet;
			// This is the extra contribution of g23:
			g23toexx = ( H.at(4,0)*( H.at(2,2)*H.at(3,4) - H.at(2,4)*H.at(3,2)) - H.at(3,0)*( H.at(2,2)*H.at(4,4) - H.at(2,4)*H.at(4,2)) + H.at(2,0)*( H.at(3,2)*H.at(4,4) - H.at(3,4)*H.at(4,2)))/subdet;
			g23toeyy = ( H.at(4,1)*( H.at(2,2)*H.at(3,4) - H.at(2,4)*H.at(3,2)) - H.at(3,1)*( H.at(2,2)*H.at(4,4) - H.at(2,4)*H.at(4,2)) + H.at(2,1)*( H.at(3,2)*H.at(4,4) - H.at(3,4)*H.at(4,2)))/subdet;
			g23tog12 = ( H.at(4,5)*( H.at(2,2)*H.at(3,4) - H.at(2,4)*H.at(3,2)) - H.at(3,5)*( H.at(2,2)*H.at(4,4) - H.at(2,4)*H.at(4,2)) + H.at(2,5)*( H.at(3,2)*H.at(4,4) - H.at(3,4)*H.at(4,2)))/subdet;
			// This is the extra contribution of g13: 
			g13toexx = ( H.at(3,0)*( H.at(2,2)*H.at(4,3) - H.at(2,3)*H.at(4,2)) - H.at(4,0)*( H.at(2,2)*H.at(3,3) - H.at(2,3)*H.at(3,2)) - H.at(2,0)*( H.at(3,2)*H.at(4,3) - H.at(3,3)*H.at(4,2)))/subdet;
			g13toeyy = ( H.at(3,1)*( H.at(2,2)*H.at(4,3) - H.at(2,3)*H.at(4,2)) - H.at(4,1)*( H.at(2,2)*H.at(3,3) - H.at(2,3)*H.at(3,2)) - H.at(2,1)*( H.at(3,2)*H.at(4,3) - H.at(3,3)*H.at(4,2)))/subdet;
			g13tog12 = ( H.at(3,5)*( H.at(2,2)*H.at(4,3) - H.at(2,3)*H.at(4,2)) - H.at(4,5)*( H.at(2,2)*H.at(3,3) - H.at(2,3)*H.at(3,2)) - H.at(2,5)*( H.at(3,2)*H.at(4,3) - H.at(3,3)*H.at(4,2)))/subdet;

			ezztoexx.reuseit(); ezztoeyy.reuseit(); ezztog12.reuseit(); g23toexx.reuseit(); g23toeyy.reuseit(); g23tog12.reuseit(); g13toexx.reuseit(); g13toeyy.reuseit(); g13tog12.reuseit();

			H = expression(3,3,{  
			H.at(0,0) + H.at(0,2)*ezztoexx+H.at(0,3)*g23toexx+H.at(0,4)*g13toexx,
			H.at(0,1) + H.at(0,2)*ezztoeyy+H.at(0,3)*g23toeyy+H.at(0,4)*g13toeyy,
			H.at(0,5) + H.at(0,2)*ezztog12+H.at(0,3)*g23tog12+H.at(0,4)*g13tog12,

			H.at(1,0) + H.at(1,2)*ezztoexx+H.at(1,3)*g23toexx+H.at(1,4)*g13toexx,
			H.at(1,1) + H.at(1,2)*ezztoeyy+H.at(1,3)*g23toeyy+H.at(1,4)*g13toeyy,
			H.at(1,5) + H.at(1,2)*ezztog12+H.at(1,3)*g23tog12+H.at(1,4)*g13tog12,

			H.at(5,0) + H.at(5,2)*ezztoexx+H.at(5,3)*g23toexx+H.at(5,4)*g13toexx,
			H.at(5,1) + H.at(5,2)*ezztoeyy+H.at(5,3)*g23toeyy+H.at(5,4)*g13toeyy,
			H.at(5,5) + H.at(5,2)*ezztog12+H.at(5,3)*g23tog12+H.at(5,4)*g13tog12
			});
		}
		
		if (myoption == "planestrain" || myoption == "planestress")
		{
			H.reuseit();
			
			expression gradu = grad(u);
			gradu.reuseit();

			expression Ei = greenlagrangestrain(gradu);
			Ei.reuseit();

			expression Si = H*Ei; 
			Si.reuseit();

			expression graddofdu = grad(dofu)-gradu;
			expression gradtfdu = grad(tfu);

			expression ei = 0.5*( graddofdu + transpose(graddofdu) + transpose(gradu)*graddofdu + transpose(graddofdu)*gradu );
			ei = expression(3,1, { entry(0,0,ei),entry(1,1,ei),2*entry(0,1,ei) });

			expression deltae = 0.5*( gradtfdu + transpose(gradtfdu) + transpose(gradu)*gradtfdu + transpose(gradtfdu)*gradu );
			deltae = expression(3,1, { entry(0,0,deltae),entry(1,1,deltae),2*entry(0,1,deltae) });

			expression deltaeta = 0.5*( transpose(graddofdu) * gradtfdu + transpose(gradtfdu)*graddofdu );
			deltaeta = expression(3,1, { entry(0,0,deltaeta),entry(1,1,deltaeta),2*entry(0,1,deltaeta) });

			prestress.reuseit();

			if (prestress.iszero())
				return -(H*ei)*deltae - Si*deltaeta - Si*deltae;
			else
				return -(H*ei)*deltae - (Si+prestress)*deltaeta - (Si+prestress)*deltae;
		}

		// If the option is not valid:
		std::cout << "Error in 'mathop' namespace: invalid option or no option provided for the 2D problem in 'predefinedelasticity'" << std::endl;
		std::cout << "Available choices are: 'planestrain', 'planestress'" << std::endl;
		abort();
	}
    if (dofu.countrows() == 3)
    {
		if (prestress.iszero() == false && (prestress.countcolumns() != 1 || prestress.countrows() != 6))
		{
		    std::cout << "Error in 'mathop' namespace: expected a 6x1 sized prestress vector (Voigt form) in 'predefinedelasticity' (set scalar 0.0 if no prestress)" << std::endl;
		    abort();
		}

        if (myoption.length() > 0)
        {
            std::cout << "Error in 'mathop' namespace: for a 3D problem the last string argument must be empty in 'predefinedelasticity'" << std::endl;
            abort();
        }

		H.reuseit();

		expression gradu = grad(u);
		gradu.reuseit();

		expression Ei = greenlagrangestrain(gradu);
		Ei.reuseit();

		expression Si = H*Ei; 
		Si.reuseit();

		expression graddofdu = grad(dofu)-gradu;
		expression gradtfdu = grad(tfu);

		expression ei = 0.5*( graddofdu + transpose(graddofdu) + transpose(gradu)*graddofdu + transpose(graddofdu)*gradu );
		ei = expression(6,1, { entry(0,0,ei),entry(1,1,ei),entry(2,2,ei),2*entry(1,2,ei),2*entry(0,2,ei),2*entry(0,1,ei) });

		expression deltae = 0.5*( gradtfdu + transpose(gradtfdu) + transpose(gradu)*gradtfdu + transpose(gradtfdu)*gradu );
		deltae = expression(6,1, { entry(0,0,deltae),entry(1,1,deltae),entry(2,2,deltae),2*entry(1,2,deltae),2*entry(0,2,deltae),2*entry(0,1,deltae) });

		expression deltaeta = 0.5*( transpose(graddofdu) * gradtfdu + transpose(gradtfdu)*graddofdu );
		deltaeta = expression(6,1, { entry(0,0,deltaeta),entry(1,1,deltaeta),entry(2,2,deltaeta),2*entry(1,2,deltaeta),2*entry(0,2,deltaeta),2*entry(0,1,deltaeta) });

		prestress.reuseit();

		if (prestress.iszero())
			return -(H*ei)*deltae - Si*deltaeta - Si*deltae;
		else
			return -(H*ei)*deltae - (Si+prestress)*deltaeta - (Si+prestress)*deltae;
    }
}

expression mathop::predefinedelectrostaticforce(expression tfu, expression E, expression epsilon)
{
    if (tfu.countcolumns() != 1)
    {
        std::cout << "Error in 'mathop' namespace: the force formula expected a column vector expression as first argument" << std::endl;
        abort();
    }

	std::vector<expression> spacederivatives(tfu.countrows());
	for (int i = 0; i < tfu.countrows(); i++)
		spacederivatives[i] = grad(tfu.at(i,0));
    
    return predefinedelectrostaticforce(spacederivatives, E, epsilon);
}

expression mathop::predefinedelectrostaticforce(std::vector<expression> dxyztfu, expression E, expression epsilon)
{
    E.reuseit();
    epsilon.reuseit();
    
    std::vector<std::vector<expression>> exprs(dxyztfu.size());
    for (int i = 0; i < dxyztfu.size(); i++)
    	exprs[i] = {dxyztfu[i]};
    
    // Scalar gradient here:
    expression gradtfu(exprs);
    
    if (gradtfu.countcolumns() == 1)
    {
        std::cout << "Error in 'mathop' namespace: the force formula is undefined for 1D displacements" << std::endl;
        abort();
    }
    
    if (gradtfu.countcolumns() == 2)
        return -( epsilon*0.5 * (pow(compx(E),2) * entry(0,0,gradtfu) - pow(compy(E),2) * entry(0,0,gradtfu) + 2 * compx(E) * compy(E) * entry(1,0,gradtfu))      +epsilon*0.5 * (-pow(compx(E),2) * entry(1,1,gradtfu) + pow(compy(E),2) * entry(1,1,gradtfu) + 2 * compy(E) * compx(E) * entry(0,1,gradtfu)) );
    if (gradtfu.countcolumns() == 3)
        return -( epsilon*0.5 * (pow(compx(E),2) * entry(0,0,gradtfu) - pow(compy(E),2) * entry(0,0,gradtfu) - pow(compz(E),2) * entry(0,0,gradtfu) + 2 * compx(E) * compy(E) * entry(1,0,gradtfu) + 2 * compx(E) * compz(E) * entry(2,0,gradtfu))      +epsilon*0.5 * (-pow(compx(E),2) * entry(1,1,gradtfu) + pow(compy(E),2) * entry(1,1,gradtfu) - pow(compz(E),2) * entry(1,1,gradtfu) + 2 * compy(E) * compx(E) * entry(0,1,gradtfu) + 2 * compy(E) * compz(E) * entry(2,1,gradtfu))      +epsilon*0.5 * (-pow(compx(E),2) * entry(2,2,gradtfu) - pow(compy(E),2) * entry(2,2,gradtfu) + pow(compz(E),2) * entry(2,2,gradtfu) + 2 * compz(E) * compx(E) * entry(0,2,gradtfu) + 2 * compz(E) * compy(E) * entry(1,2,gradtfu)) );
}

expression mathop::predefinedmagnetostaticforce(expression tfu, expression H, expression mu)
{
	return predefinedelectrostaticforce(tfu, H, mu);
}

expression mathop::predefinedmagnetostaticforce(std::vector<expression> dxyztfu, expression H, expression mu)
{
	return predefinedelectrostaticforce(dxyztfu, H, mu);
}

expression mathop::predefinedstokes(expression dofv, expression tfv, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant, bool isviscosityconstant)
{
    int problemdimension = universe::mymesh->getmeshdimension();

    if (problemdimension < 2)
    {
        std::cout << "Error in 'mathop' namespace: 'predefinedstokes' is only allowed on 2D and 3D geometries" << std::endl;
        abort();
    }
    if (universe::isaxisymmetric)
        problemdimension++;
    if (dofv.countcolumns() != 1 || dofv.countrows() != problemdimension || tfv.countcolumns() != 1 || tfv.countrows() != problemdimension || not(dofp.isscalar()) || not(tfp.isscalar()) || not(mu.isscalar()) || not(rho.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinedstokes'" << std::endl;
        abort();
    }

    expression output = predefinedmassconservation(dofv, tfp, rho, dtrho, gradrho, includetimederivs, isdensityconstant);
    
    if (includetimederivs)
        output = output - rho*dt(dofv)*tfv;

    return ( output - grad(dofp)*tfv + predefinedviscousforce(dofv, tfv, mu, isdensityconstant, isviscosityconstant) );
}

expression mathop::predefinednavierstokes(expression dofv, expression tfv, expression v, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant, bool isviscosityconstant)
{
    int problemdimension = universe::mymesh->getmeshdimension();

    if (problemdimension < 2)
    {
        std::cout << "Error in 'mathop' namespace: 'predefinednavierstokes' is only allowed on 2D and 3D geometries" << std::endl;
        abort();
    }
    if (universe::isaxisymmetric)
        problemdimension++;
    if (dofv.countcolumns() != 1 || dofv.countrows() != problemdimension || tfv.countcolumns() != 1 || tfv.countrows() != problemdimension || v.countcolumns() != 1 || v.countrows() != problemdimension || not(dofp.isscalar()) || not(tfp.isscalar()) || not(mu.isscalar()) || not(rho.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinednavierstokes'" << std::endl;
        abort();
    }
    
    expression output = predefinedmassconservation(dofv, tfp, rho, dtrho, gradrho, includetimederivs, isdensityconstant);
    
    if (includetimederivs)
        output = output - rho*dt(dofv)*tfv;
        
    return ( output - grad(dofp)*tfv + predefinedviscousforce(dofv, tfv, mu, isdensityconstant, isviscosityconstant) +  predefinedinertialforce(dofv, tfv, v, rho) );
}




