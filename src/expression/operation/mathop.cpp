#include "mathop.h"


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

    // Make sure the filename includes the extension:
    if (filename.size() < 5 || filename.substr(filename.size()-4,4) != ".pos")
    {
        std::cout << "Error in 'mathop' namespace: cannot write to file '" << filename << "' (unknown or missing file extension)" << std::endl;
        abort();
    }
    // Remove the extension:
    filename = filename.substr(0, filename.size()-4);

    // Write the header:
    gmshinterface::openview(filename + ".pos", filename, 0, true);
    // Write the data:
    if (isscalar)
        gmshinterface::appendtoview(filename + ".pos", 0, densematrix(n,1,xcoords), densematrix(n,1,ycoords), densematrix(n,1,zcoords), densematrix(n,1,compxevals));
    else
        gmshinterface::appendtoview(filename + ".pos", 0, densematrix(n,1,xcoords), densematrix(n,1,ycoords), densematrix(n,1,zcoords), densematrix(n,1,compxevals), densematrix(n,1,compyevals), densematrix(n,1,compzevals));
    // Close view:
    gmshinterface::closeview(filename + ".pos");
}
    
    
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
        std::cout << "Error in 'mathop' namespace: can only take the gradient of scalars and column vectors" << std::endl;
        abort();
    }

    int problemdimension = universe::mymesh->getmeshdimension();

    std::vector<expression> myexprs = {};
    for (int i = 0; i < problemdimension; i++)
    {
        for (int comp = 0; comp < input.countrows(); comp++)
            myexprs.push_back(input.spacederivative(i+1).at(comp,0));
    }
    
    return expression(problemdimension, input.countrows(), myexprs);
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
        switch (input.countrows())
        {
            case 1:
                return expression(3,1,{0, 0, 0});
            case 2:
                return expression(3,1,{0, 0, dx(compy(input))-dy(compx(input))});
            case 3:
                return expression(3,1,{dy(compz(input))-dz(compy(input)), dz(compx(input))-dx(compz(input)), dx(compy(input))-dy(compx(input))});
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
    return expression(1,1, {term11});
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
 
vec mathop::solve(mat A, vec b)
{
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



////////// PREDEFINED OPERATORS

expression mathop::strain(expression input)
{
    if ((input.countrows() != 2 && input.countrows() != 3) || input.countcolumns() != 1)
    {
        std::cout << "Error in 'mathop' namespace: can only compute the strains of a 2x1 or 3x1 column vector" << std::endl;
        abort();
    }
	if (input.countrows() == 2)
		return expression(3,1,{dx(compx(input)), dy(compy(input)), dy(compx(input)) + dx(compy(input))});
	if (input.countrows() == 3)
		return expression(6,1,{dx(compx(input)), dy(compy(input)), dz(compz(input)), dz(compy(input)) + dy(compz(input)), dx(compz(input)) + dz(compx(input)), dx(compy(input)) + dy(compx(input))});
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
        std::cout << "Error in 'mathop' namespace: first two arguments in 'predefinedelasticity' must be both either 2x1 or 3x1 vectors" << std::endl;
        abort();
    }
    if (dofu.countrows() == 2)
	{
		if (myoption == "planestrain")
		{
			expression Hplanestrain = expression(3,3, {H.at(0,0),H.at(0,1),H.at(0,5), H.at(1,0),H.at(1,1),H.at(1,5), H.at(5,0),H.at(5,1),H.at(5,5) });
			return -transpose( Hplanestrain *strain(dofu) )*strain(tfu);
		}
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

			expression Hplanestress(3,3,{  
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

			return -transpose( Hplanestress *strain(dofu) )*strain(tfu);
		}
		
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
        return -transpose( H*strain(dofu) ) * strain(tfu);
    }
}

expression mathop::predefinedelectrostaticforce(expression gradtfu, expression gradv, expression epsilon)
{
    expression E = gradv;
    E.reuseit();
    
    if (gradtfu.countcolumns() == 1)
    {
        std::cout << "Error in 'mathop' namespace: 'predefinedelectrostaticforce' is undefined for 1D displacements" << std::endl;
        abort();
    }
    if (gradtfu.countcolumns() == 2)
        return -( epsilon*0.5 * (pow(compx(E),2) * entry(0,0,gradtfu) - pow(compy(E),2) * entry(0,0,gradtfu) + 2 * compx(E) * compy(E) * entry(1,0,gradtfu))      +epsilon*0.5 * (-pow(compx(E),2) * entry(1,1,gradtfu) + pow(compy(E),2) * entry(1,1,gradtfu) + 2 * compy(E) * compx(E) * entry(0,1,gradtfu)) );
    if (gradtfu.countcolumns() == 3)
        return -( epsilon*0.5 * (pow(compx(E),2) * entry(0,0,gradtfu) - pow(compy(E),2) * entry(0,0,gradtfu) - pow(compz(E),2) * entry(0,0,gradtfu) + 2 * compx(E) * compy(E) * entry(1,0,gradtfu) + 2 * compx(E) * compz(E) * entry(2,0,gradtfu))      +epsilon*0.5 * (-pow(compx(E),2) * entry(1,1,gradtfu) + pow(compy(E),2) * entry(1,1,gradtfu) - pow(compz(E),2) * entry(1,1,gradtfu) + 2 * compy(E) * compx(E) * entry(0,1,gradtfu) + 2 * compy(E) * compz(E) * entry(2,1,gradtfu))      +epsilon*0.5 * (-pow(compx(E),2) * entry(2,2,gradtfu) - pow(compy(E),2) * entry(2,2,gradtfu) + pow(compz(E),2) * entry(2,2,gradtfu) + 2 * compz(E) * compx(E) * entry(0,2,gradtfu) + 2 * compz(E) * compy(E) * entry(1,2,gradtfu)) );
}

