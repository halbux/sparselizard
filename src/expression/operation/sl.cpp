#include "sl.h"
#include <random>
#include "rawpoint.h"
#include "rawline.h"
#include "rawsurface.h"
#include "rawvolume.h"
#include "slmpi.h"


int sl::getversion(void)
{
    return 202111;
}

int sl::getsubversion(void)
{
    return 10;
}

std::string sl::getversionname(void)
{
    return "elegant elk";
}

void sl::printversion(void)
{
    std::cout << "sparselizard " << std::to_string(getversion()) << "." << std::to_string(getsubversion()) << " ('" << getversionname() << "')." << std::endl;
}

void sl::setmaxnumthreads(int mnt)
{
    if (mnt <= 0)
    {
        std::cout << "Error in 'sl' namespace: cannot set a negative or zero max num threads" << std::endl;
        abort();
    }

    universe::setmaxnumthreads(mnt);
}

int sl::getmaxnumthreads(void)
{
    return universe::getmaxnumthreads();
}

double sl::getpi(void)
{
    return 3.1415926535897932384;
}

double sl::getrandom(void)
{
    // Will be used to obtain a seed for the random number engine:
    std::random_device rd;
    // Standard mersenne twister engine seeded with rd():
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    return dis(gen);
}

int sl::selectunion(std::vector<int> physregs)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined(physregs);
    
    universe::getrawmesh()->getoriginalmeshpointer()->getphysicalregions()->createunion(physregs, false);
    return (universe::getrawmesh()->getphysicalregions())->createunion(physregs, false);
}

int sl::selectintersection(std::vector<int> physregs)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined(physregs);
    
    universe::getrawmesh()->getoriginalmeshpointer()->getphysicalregions()->createintersection(physregs, false);
    return (universe::getrawmesh()->getphysicalregions())->createintersection(physregs, false);
}

int sl::selectall(void)
{
    universe::getrawmesh()->getoriginalmeshpointer()->getphysicalregions()->createunionofall(false);
    return (universe::getrawmesh()->getphysicalregions())->createunionofall(false);
}

bool sl::isdefined(int physreg)
{
    return (universe::getrawmesh()->getphysicalregions()->getindex(physreg) != -1);
}

bool sl::isempty(int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});

    std::vector<bool> defin = universe::getrawmesh()->getphysicalregions()->get(physreg)->getdefinition();
    
    for (int i = 0; i < defin.size(); i++)
    {
        if (defin[i] == true)
            return false;
    }
    
    return true;
}

bool sl::isinside(int physregtocheck, int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physregtocheck,physreg});
    
    std::vector<bool> tocheckdefin = universe::getrawmesh()->getphysicalregions()->get(physregtocheck)->getdefinition();
    std::vector<bool> defin = universe::getrawmesh()->getphysicalregions()->get(physreg)->getdefinition();
    
    // All disjoint regions must be included:
    for (int i = 0; i < tocheckdefin.size(); i++)
    {
        if (tocheckdefin[i] == true && defin[i] == false)
            return false;
    }

    return true;
}

bool sl::istouching(int physregtocheck, int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physregtocheck,physreg});
    
    std::vector<bool> tocheckdefin = universe::getrawmesh()->getphysicalregions()->get(physregtocheck)->getdefinition();
    std::vector<bool> defin = universe::getrawmesh()->getphysicalregions()->get(physreg)->getdefinition();
    
    // At least one disjoint region must be included:
    for (int i = 0; i < tocheckdefin.size(); i++)
    {
        if (tocheckdefin[i] == true && defin[i] == true)
            return true;
    }

    return false;
}

void sl::printvector(std::vector<double> input)
{
    std::cout << "Vector size is " << input.size() << std::endl;
    for (int i = 0; i < input.size(); i++)
        std::cout << input[i] << " ";
    std::cout << std::endl;
}

void sl::printvector(std::vector<int> input)
{
    std::cout << "Vector size is " << input.size() << std::endl;
    for (int i = 0; i < input.size(); i++)
        std::cout << input[i] << " ";
    std::cout << std::endl;
}

void sl::printvector(std::vector<bool> input)
{
    std::cout << "Vector size is " << input.size() << std::endl;
    for (int i = 0; i < input.size(); i++)
        std::cout << input[i] << " ";
    std::cout << std::endl;
}

void sl::writevector(std::string filename, std::vector<double> towrite, char delimiter, bool writesize)
{
    if (towrite.size() == 0)
        return;

    // 'file' cannot take a std::string argument --> filename.c_str():
    std::ofstream name (filename.c_str());
    if (name.is_open())
    {
        if (writesize)
            name << towrite.size() << delimiter;

        // To write all doubles with enough digits to the file:
        name << std::setprecision(17);

        for (int i = 0; i < towrite.size()-1; i++)
            name << towrite[i] << delimiter;
        name << towrite[towrite.size()-1];

        name.close();
    }
    else
    {
        std::cout << "Unable to write vector to file " << filename << " or file not found" << std::endl;
        abort();
    }
}

std::vector<double> sl::loadvector(std::string filename, char delimiter, bool sizeincluded)
{
    std::vector<double> output = {};

    std::string currentline;

    // 'file' cannot take a std::string argument --> filename.c_str():
    std::ifstream name (filename.c_str());
    if (name.is_open())
    {
        if (sizeincluded)
        {
            std::getline(name, currentline, delimiter);
            gentools::osclean(currentline);
            output.resize(std::stoi(currentline));
        }

        int index = 0;
        while (std::getline(name, currentline, delimiter))
        {
            gentools::osclean(currentline);
            if (sizeincluded)
                output[index] = std::stod(currentline);
            else
                output.push_back(std::stod(currentline));

            index++;
        }

        name.close();
    }
    else
    {
        std::cout << "Unable to load vector from file " << filename << " or file not found" << std::endl;
        abort();
    }

    return output;
}

#ifndef HAVE_GMSH
std::string sl::allpartition(std::string meshfile)
{
    if (slmpi::count() == 1)
        return meshfile;
    
    std::cout << "Error in 'sl' namespace: GMSH API is required to partition the mesh" << std::endl;
    abort();
}
#endif
#ifdef HAVE_GMSH
#include "gmsh.h"
std::string sl::allpartition(std::string meshfile)
{
    int rank = slmpi::getrank();
    int numranks = slmpi::count();
    
    if (numranks == 1)
        return meshfile;
    
    if (rank == 0)
    {
        gmsh::initialize();
        gmsh::open(meshfile);
        // Unfortunately dropping the global info and saving in format 2 is not allowed
        // gmsh::option::setNumber("Mesh.Format", 2);
        gmsh::option::setNumber("Mesh.PartitionSplitMeshFiles", 1);
        gmsh::option::setNumber("Mesh.PartitionCreateTopology", 0);
        gmsh::model::mesh::partition(numranks);
        gmsh::write(meshfile);
        gmsh::finalize();
    }
    // Wait for rank 0 to finish:
    slmpi::barrier();
    
    return (meshfile.substr(0, meshfile.size()-4)+"_"+std::to_string(rank+1)+".msh");
}
#endif

expression sl::norm(expression expr)
{
    if (expr.isscalar())
        return abs(expr);

    expression mynorm;
    
    for (int i = 0; i < expr.countrows(); i++)
    {
        for (int j = 0; j < expr.countcolumns(); j++)
        {
            if (i == 0 && j == 0)
                mynorm = pow(expr.at(i,j),2);
            else
                mynorm = mynorm + pow(expr.at(i,j),2);
        }
    }

    return sqrt(mynorm);
}

expression sl::normal(void)
{
    return getnormal(-1);
}

expression sl::normal(int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return getnormal(physreg);
}

expression sl::getnormal(int physreg)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();

    if (physreg >= 0)
    {
        int elementdimension = universe::getrawmesh()->getphysicalregions()->get(physreg)->getelementdimension();
        if (elementdimension >= 0 && elementdimension != problemdimension)
        {
            std::cout << "Error in 'sl' namespace: normal cannot point outward of the " << elementdimension << "D region provided (should be " << problemdimension << "D)" << std::endl;
            abort();
        }
    }
    
    expression output;

    expression expr;
    if (problemdimension == 1)
        output = 1.0;
    if (problemdimension == 2)
    {
        expression mynorm = sqrt(expr.invjac(0,1)*expr.invjac(0,1)+expr.invjac(1,1)*expr.invjac(1,1));
        mynorm.reuseit();
        if (universe::isaxisymmetric)
            output = array3x1(expr.invjac(0,1), expr.invjac(1,1), 0)/mynorm;
        else
            output = array2x1(expr.invjac(0,1), expr.invjac(1,1))/mynorm;
    }
    if (problemdimension == 3)
    {
        expression mynorm = sqrt(expr.invjac(0,2)*expr.invjac(0,2)+expr.invjac(1,2)*expr.invjac(1,2)+expr.invjac(2,2)*expr.invjac(2,2));
        mynorm.reuseit();
        output = array3x1(expr.invjac(0,2), expr.invjac(1,2), expr.invjac(2,2))/mynorm;
    }
    
    if (physreg >= 0)
    {
        std::shared_ptr<oporientation> op(new oporientation(physreg));
        expression exprorient(op);
        exprorient.reuseit();
        output = exprorient*output;
    }

    output.reuseit();

    return output;
}

expression sl::tangent(void)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();

    expression expr;
    if (problemdimension == 1)
        return 1.0;
    if (problemdimension == 2)
    {
        expression mynorm = sqrt(expr.jac(0,0)*expr.jac(0,0)+expr.jac(0,1)*expr.jac(0,1));
        mynorm.reuseit();
        if (universe::isaxisymmetric)
            return array3x1(expr.jac(0,0), expr.jac(0,1), 0)/mynorm;
        else
            return array2x1(expr.jac(0,0), expr.jac(0,1))/mynorm;
    }
    if (problemdimension == 3)
    {
        expression mynorm = sqrt(expr.jac(0,0)*expr.jac(0,0)+expr.jac(0,1)*expr.jac(0,1)+expr.jac(0,2)*expr.jac(0,2));
        mynorm.reuseit();
        return array3x1(expr.jac(0,0), expr.jac(0,1), expr.jac(0,2))/mynorm;
    }
    
    abort(); // fix return warning
}

void sl::scatterwrite(std::string filename, std::vector<double> xcoords, std::vector<double> ycoords, std::vector<double> zcoords, std::vector<double> compxevals, std::vector<double> compyevals, std::vector<double> compzevals)
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
        std::cout << "Error in 'sl' namespace: size of 'scatterwrite' arguments do not match" << std::endl;
        abort();
    }

    iodata datatowrite(1, 1, isscalar, {});
    datatowrite.addcoordinates(0, densemat(n,1,xcoords), densemat(n,1,ycoords), densemat(n,1,zcoords));
    if (isscalar)
        datatowrite.adddata(0, {densemat(n,1,compxevals)});
    else
        datatowrite.adddata(0, {densemat(n,1,compxevals), densemat(n,1,compyevals), densemat(n,1,compzevals)});

    iointerface::writetofile(filename, datatowrite);
}


void sl::setaxisymmetry(void)
{
    // Make sure the call is done before loading the mesh:
    if (universe::myrawmesh != NULL)
    {
        std::cout << "Error in 'sl' namespace: 'setaxisymmetry' must be called before loading the mesh" << std::endl;
        abort();
    }
    universe::isaxisymmetric = true;
}

void sl::setfundamentalfrequency(double f) { universe::fundamentalfrequency = f; }
void sl::settime(double t) { universe::currenttimestep = t; }
double sl::gettime(void) { return universe::currenttimestep; }

expression sl::meshsize(int integrationorder)
{
    std::shared_ptr<opmeshsize> op(new opmeshsize(integrationorder));
    return expression(op);
}

expression sl::fieldorder(field input, double alpha, double absthres)
{
    std::shared_ptr<rawfield> rf = input.getpointer();
    
    if (rf->gettypename() != "h1" && rf->gettypename() != "hcurl" && rf->gettypename() != "h1d0" && rf->gettypename() != "h1d1" && rf->gettypename() != "h1d2" && rf->gettypename() != "h1d3")
    {
        std::cout << "Error in 'sl' namespace: field provided to 'fieldorder' must be of type 'h1', 'h1d' or 'hcurl' (was '" << rf->gettypename() << "')" << std::endl;
        abort();
    }
    
    std::shared_ptr<opfieldorder> op(new opfieldorder(rf->getsons(), alpha, absthres));

    return expression(op);
}

expression sl::getharmonic(int harmnum, expression input, int numfftharms)
{
    if (harmnum <= 0)
    {
        std::cout << "Error in 'sl' namespace: in 'getharmonic' cannot have a negative or zero harmonic" << std::endl;
        abort();
    }
    
    if (input.iszero())
        return input;
    
    return moveharmonic({harmnum}, {1}, input, numfftharms);
}

expression sl::makeharmonic(std::vector<int> harms, std::vector<expression> exprs)
{
    if (harms.size() == 0)
    {
        std::cout << "Error in 'sl' namespace: in 'makeharmonic' expected at least one harmonic as argument" << std::endl;
        abort();
    }
    if (harms.size() != exprs.size())
    {
        std::cout << "Error in 'sl' namespace: in 'makeharmonic' the number of harmonics and expressions do not match" << std::endl;
        abort();
    }
    
    int m = exprs[0].countrows();
    int n = exprs[0].countcolumns();
    
    std::vector<int> alldisjregs(universe::getrawmesh()->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    
    for (int i = 0; i < harms.size(); i++)
    {
        if (harms[i] <= 0)
        {
            std::cout << "Error in 'sl' namespace: in 'makeharmonic' cannot have a negative or zero harmonic" << std::endl;
            abort();
        }
        if (exprs[i].countrows() != m || exprs[i].countcolumns() != n)
        {
            std::cout << "Error in 'sl' namespace: in 'makeharmonic' all expressions should have the same dimension" << std::endl;
            abort();
        }
        if (not(exprs[i].isharmonicone(alldisjregs)))
        {
            std::cout << "Error in 'sl' namespace: in 'makeharmonic' cannot have multiharmonic expressions as argument (only constant harmonic 1)" << std::endl;
            abort();
        }
    }
    
    if (harms.size() == 1)
        return moveharmonic({1}, {harms[0]}, exprs[0]);
    
    std::vector<expression> outexprs(m*n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::vector<std::shared_ptr<operation>> ops(harms.size());
            for (int k = 0; k < harms.size(); k++)
                ops[k] = moveharmonic({1},{harms[k]}, exprs[k].at(i,j)).getoperationinarray(0,0);
            std::shared_ptr<opsum> op(new opsum(ops));
            outexprs[i*n+j] = expression(op);
        }
    }
    
    return expression(m, n, outexprs);
}

expression sl::moveharmonic(std::vector<int> origharms, std::vector<int> destharms, expression input, int numfftharms)
{
    if (origharms.size() == 0)
    {
        std::cout << "Error in 'sl' namespace: in 'moveharmonic' expected at least one harmonic as argument" << std::endl;
        abort();
    }
    if (origharms.size() != destharms.size())
    {
        std::cout << "Error in 'sl' namespace: in 'moveharmonic' the number of origin and destination harmonics do not match" << std::endl;
        abort();
    }
    for (int i = 0; i < origharms.size(); i++)
    {
        if (origharms[i] <= 0 || destharms[i] <= 0)
        {
            std::cout << "Error in 'sl' namespace: in 'moveharmonic' cannot have a negative or zero harmonic" << std::endl;
            abort();
        }
    }
    
    int m = input.countrows();
    int n = input.countcolumns();
    
    std::vector<expression> exprs(m*n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::shared_ptr<operation> opin = input.getoperationinarray(i,j);
            
            if (opin->isdofincluded() || opin->istfincluded())
            {
                std::cout << "Error in 'sl' namespace: in 'moveharmonic' expected an argument expression without dof or tf" << std::endl;
                abort();
            }

            std::shared_ptr<opharmonic> op(new opharmonic(origharms, destharms, opin, numfftharms));
            exprs[i*n+j] = expression(op);
        }
    }

    return expression(m,n,exprs);
}

std::vector<double> sl::gettotalforce(int physreg, expression* meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    std::vector<int> alldisjregs(universe::getrawmesh()->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(EorH.isharmonicone(alldisjregs)) || not(epsilonormu.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'sl' namespace: cannot have a multiharmonic argument in the total force calculation" << std::endl;
        abort();
    }
    
    int wholedomain = universe::getrawmesh()->getphysicalregions()->createunionofall();

    int numcomps = universe::getrawmesh()->getmeshdimension();
    if (universe::isaxisymmetric)
        numcomps++;
    if (numcomps == 1)
    {
        std::cout << "Error in 'sl' namespace: force calculation formula is undefined in 1D" << std::endl;
        abort();
    }
        
    std::vector<std::string> tn = {"","h1","h1xy","h1xyz"};

    field f(tn[numcomps]);
    f.setorder(wholedomain, 1); // must be order 1
    
    formulation forcecalc;
    if (meshdeform == NULL)
        forcecalc += integral(wholedomain, -predefinedelectrostaticforce(tf(f, physreg), EorH, epsilonormu), extraintegrationorder);
    else
        forcecalc += integral(wholedomain, *meshdeform, -predefinedelectrostaticforce(tf(f, physreg), EorH, epsilonormu), extraintegrationorder);
    forcecalc.generate();
    vec fv = forcecalc.b();
    
    // First x then y then z component in structure:
    std::vector<double> output(numcomps);
    
    int siz = fv.size()/numcomps;
    for (int c = 0; c < numcomps; c++)
    {
        // totalforce = integral(fdensity*1) = integral(fdensity*sumi(Ni)) = sumi(integral(fdensity*Ni))
        // where Ni includes all 'h1' nodal shape functions (their sum equals one). 
        densemat vecvals = fv.getvalues(indexmat(siz, 1, c*siz, 1));
        output[c] = vecvals.sum();
    }
    
    if (universe::isaxisymmetric)
        output = {0, 2.0*getpi()*output[1], 0};
    
    universe::getrawmesh()->getphysicalregions()->remove({wholedomain}, false);
    
    return output;
}

std::vector<double> sl::gettotalforce(int physreg, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return gettotalforce(physreg, NULL, EorH, epsilonormu, extraintegrationorder);
}

std::vector<double> sl::gettotalforce(int physreg, expression meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return gettotalforce(physreg, &meshdeform, EorH, epsilonormu, extraintegrationorder);
}

std::vector<double> sl::printtotalforce(int physreg, expression* meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    std::vector<double> totforce = gettotalforce(physreg, meshdeform, EorH, epsilonormu, extraintegrationorder);

    std::cout << "Total force on region " << physreg << " is ";
    std::vector<std::string> compstr = {"x","y","z"};

    for (int c = 0; c < totforce.size(); c++)
    {
        std::cout << "f" << compstr[c] << " = " << totforce[c];
        if (c != totforce.size()-1)
            std::cout << ", ";
    }
    
    std::vector<std::string> unitstr = {"",""," N per unit depth"," N"};
    std::cout << unitstr[totforce.size()] << std::endl;
    
    return totforce;
}

std::vector<double> sl::printtotalforce(int physreg, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return printtotalforce(physreg, NULL, EorH, epsilonormu, extraintegrationorder);
}

std::vector<double> sl::printtotalforce(int physreg, expression meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return printtotalforce(physreg, &meshdeform, EorH, epsilonormu, extraintegrationorder);
}

void sl::setphysicalregionshift(int shiftamount) { universe::physregshift = shiftamount; }

void sl::writeshapefunctions(std::string filename, std::string sftypename, int elementtypenumber, int maxorder, bool allorientations)
{
    if (elementtypenumber == 7)
    {
        std::cout << "Error in 'sl' namespace: cannot write shape functions for pyramids (non-polynomial)" << std::endl;
        abort();
    }
    if (elementtypenumber > 7)
    {
        std::cout << "Error in 'sl' namespace: element type number must be between 0 and 7" << std::endl;
        abort();
    }
    
    element myelement(elementtypenumber);
    int numnodes = myelement.countnodes();

    lagrangeformfunction slff(elementtypenumber, 1, {});
    std::vector<double> scoords = slff.getnodecoordinates();
    lagrangeformfunction lff(elementtypenumber, maxorder+1, {});
    std::vector<double> coords = lff.getnodecoordinates();

    std::shared_ptr<hierarchicalformfunction> hff = selector::select(elementtypenumber, sftypename);
    int numcomp = hff->countcomponents();
    std::vector<densemat> ffvals(numcomp);

    // Prepare the element x, y and z coordinates:
    densemat x(1,numnodes), y(1,numnodes), z(1,numnodes);
    double* xptr = x.getvalues();
    double* yptr = y.getvalues();
    double* zptr = z.getvalues();
    for (int i = 0; i < numnodes; i++)
    {
        xptr[i] = scoords[3*i+0];
        yptr[i] = scoords[3*i+1];
        zptr[i] = scoords[3*i+2];
    }
    // Evaluate all form functions:
    hierarchicalformfunctioncontainer hffc = hff->evalat(maxorder);
    hffc.evaluate(coords);
    
    std::vector<std::string> dimtype = {"node","edge","face","volume"};

    // Loop on all form functions:
    hierarchicalformfunctioniterator myiterator(sftypename, elementtypenumber, maxorder);
    for (int ff = 0; ff < myiterator.count(); ff++)
    {    
        int h = myiterator.getformfunctionorder();
        int i = myiterator.getdimension();
        int j = myiterator.getnodeedgefacevolumeindex();
        
        int l = myiterator.getformfunctionindexincurrentorderinnodeedgefacevolume();
        int ffnum = myiterator.getoverallformfunctionindex();
     
        int subtype = myiterator.getassociatedelementtype();   
        int numorientations = orientation::countorientations(subtype);
        if (allorientations == false)
            numorientations = 1;
        for (int o = 0; o < numorientations; o++)
        {
            for (int c = 0; c < numcomp; c++)
                ffvals[c] = hffc.tomatrix(h, i, j, o, l, 0, c);
                 
            iodata datatowrite(maxorder+1, 1, numcomp == 1, {});
                 
            datatowrite.addcoordinates(elementtypenumber, x, y, z);
            datatowrite.adddata(elementtypenumber, ffvals);
                
            std::string ori = "";
            if (allorientations)
                ori = "_orientation_" + std::to_string(o);
                
            iointerface::writetofile(filename, datatowrite, std::to_string(ffnum+1000) + "_" + std::to_string(l) + "th_of_order_" + std::to_string(h) + "_at_" + dimtype[i] + "_" + std::to_string(j) + ori);
        }
        
        myiterator.next();
    }    
}

expression sl::t(void) { expression exp; return exp.time(); }

void sl::grouptimesteps(std::string filename, std::vector<std::string> filestogroup, std::vector<double> timevals)
{
    iointerface::grouptimesteps(filename, filestogroup, timevals);
}

void sl::grouptimesteps(std::string filename, std::string fileprefix, int firstint, std::vector<double> timevals)
{
    int numsteps = timevals.size();

    std::vector<std::string> filestogroup(numsteps);
    for (int i = 0; i < numsteps; i++)
        filestogroup[i] = fileprefix+std::to_string(firstint+i)+".vtu";

    iointerface::grouptimesteps(filename, filestogroup, timevals);
}

std::vector<std::vector<shape>> sl::loadshape(std::string meshfile)
{
    std::shared_ptr<rawmesh> loadedmesh(new rawmesh());
    
    std::string tool, source;
    gentools::splitatcolon(meshfile, tool, source);
    if (tool.size() == 0)
        tool = "native";
        
    loadedmesh->readfromfile(tool, source);
    
    nodes* loadednodes = loadedmesh->getnodes();
    std::vector<double>* nodecoords = loadednodes->getcoordinates();
    int totalnumnodes = nodecoords->size()/3;
    elements* loadedelems = loadedmesh->getelements();
    physicalregions* loadedphysregs = loadedmesh->getphysicalregions();

    int curvatureorder = loadedelems->getcurvatureorder();
    
    std::vector<int> physregnums = loadedphysregs->getallnumbers();
    int numphysregs = physregnums.size();
        
    std::vector<std::vector<shape>> output(4,std::vector<shape>(0));
    for (int i = 0; i < numphysregs; i++)
    {
        int curphysreg = physregnums[i];
        physicalregion* pr = loadedphysregs->get(curphysreg);
        int physregdim = pr->getelementdimension();
        
        // Get the list of all elements in the physical region:
        std::vector<std::vector<int>>* curelemlist = pr->getelementlist();
        
        // Get the elements with all their nodes:
        std::vector<std::vector<int>> nodesinelements(8);
        for (int j = 0; j < 8; j++)
        {
            int numelems = curelemlist->at(j).size();
            element myelem(j, curvatureorder);
            int numnodes = myelem.countcurvednodes();
            nodesinelements[j].resize(numnodes*numelems);
            for (int k = 0; k < numelems; k++)
            {
                int curelem = curelemlist->at(j)[k];
                for (int l = 0; l < numnodes; l++)
                    nodesinelements[j][numnodes*k+l] = loadedelems->getsubelement(0,j,curelem,l);
            }
        }
        // Create a vector to renumber all node coordinates from 0 up consecutively:
        std::vector<int> renumbernodes(totalnumnodes,-1);
        int nodeindex = 0;
        for (int j = 0; j < 8; j++)
        {
            for (int k = 0; k < nodesinelements[j].size(); k++)
            {
                int curnode = nodesinelements[j][k];
                if (renumbernodes[curnode] == -1)
                {
                    renumbernodes[curnode] = nodeindex;
                    nodeindex++;
                }
                nodesinelements[j][k] = renumbernodes[curnode];
            }
        }
        // Create the vector of node coordinates:
        std::vector<double> nodecoordinates(3*nodeindex);
        for (int j = 0; j < renumbernodes.size(); j++)
        {
            if (renumbernodes[j] != -1)
            {
                nodecoordinates[3*renumbernodes[j]+0] = nodecoords->at(3*j+0);
                nodecoordinates[3*renumbernodes[j]+1] = nodecoords->at(3*j+1);
                nodecoordinates[3*renumbernodes[j]+2] = nodecoords->at(3*j+2);
            }
        }
        
        shape curshape;
        if (physregdim == 0)
            curshape = shape(std::shared_ptr<rawpoint>(new rawpoint(curphysreg, nodecoordinates, nodesinelements, curvatureorder)));
        if (physregdim == 1)
            curshape = shape(std::shared_ptr<rawline>(new rawline(curphysreg, nodecoordinates, nodesinelements, curvatureorder)));
        if (physregdim == 2)
            curshape = shape(std::shared_ptr<rawsurface>(new rawsurface(curphysreg, nodecoordinates, nodesinelements, curvatureorder)));
        if (physregdim == 3)
            curshape = shape(std::shared_ptr<rawvolume>(new rawvolume(curphysreg, nodecoordinates, nodesinelements, curvatureorder)));
    
        output[physregdim].push_back(curshape);
    }
    
    return output;
}

void sl::settimederivative(vec dtx)
{
    universe::xdtxdtdtx = {{},{dtx},{}};
}

void sl::settimederivative(vec dtx, vec dtdtx)
{
    universe::xdtxdtdtx = {{},{dtx},{dtdtx}};
}

expression sl::dx(expression input) { return input.spacederivative(1); }
expression sl::dy(expression input) { return input.spacederivative(2); }
expression sl::dz(expression input) { return input.spacederivative(3); }

expression sl::dt(expression input) { return input.timederivative(1); }
expression sl::dtdt(expression input) { return input.timederivative(2); }
expression sl::dtdtdt(expression input) { return input.timederivative(3); }
expression sl::dtdtdtdt(expression input) { return input.timederivative(4); }

expression sl::dt(expression input, double initdt, double initdtdt)
{
    return transientdtapprox(1, input, initdt, initdtdt);
}

expression sl::dtdt(expression input, double initdt, double initdtdt)
{
    return transientdtapprox(2, input, initdt, initdtdt);
}
    
expression sl::transientdtapprox(int dtorder, expression input, double initdt, double initdtdt)
{
    int m = input.countrows(), n = input.countcolumns();
    std::vector<expression> exprs(m*n);
    
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::shared_ptr<opdtapprox> op(new opdtapprox(dtorder, input.getoperationinarray(i, j), initdt, initdtdt));
            
            universe::opdtapproxes.push_back(op);
            
            exprs[i*n+j] = expression(op);
        }
    }
    
    return expression(m, n, exprs);
}

expression sl::sin(expression input) { return input.sin(); }
expression sl::cos(expression input) { return input.cos(); }
expression sl::tan(expression input) { return input.tan(); }
expression sl::asin(expression input) { return input.asin(); }
expression sl::acos(expression input) { return input.acos(); }
expression sl::atan(expression input) { return input.atan(); }
expression sl::abs(expression input) { return input.abs(); }
expression sl::sqrt(expression input) { return pow(input, 0.5); }
expression sl::log10(expression input) { return input.log10(); }
expression sl::pow(expression base, expression exponent) { return base.pow(exponent); }
expression sl::exp(expression input) { return pow(2.7182818284590452353, input); }
expression sl::mod(expression input, double modval) { return input.mod(modval); }

expression sl::ifpositive(expression condexpr, expression trueexpr, expression falseexpr)
{
    expression output(condexpr, trueexpr, falseexpr);
    return output;
}

expression sl::andpositive(std::vector<expression> exprs)
{
    if (exprs.size() == 0)
    {
        std::cout << "Error in 'sl' namespace: cannot call andpositive on an empty vector of expressions" << std::endl;
        abort();
    }

    expression output(exprs[exprs.size()-1], 1, -1);

    for (int i = exprs.size()-2; i >= 0; i--)
        output = expression(exprs[i], output, -1);

    return output;
}

expression sl::orpositive(std::vector<expression> exprs)
{
    if (exprs.size() == 0)
    {
        std::cout << "Error in 'sl' namespace: cannot call orpositive on an empty vector of expressions" << std::endl;
        abort();
    }

    expression output(exprs[exprs.size()-1], 1, -1);

    for (int i = exprs.size()-2; i >= 0; i--)
        output = expression(exprs[i], 1, output);

    return output;
}

expression sl::max(expression a, expression b)
{
    a.reuseit(); b.reuseit();
    
    return ifpositive(a-b, a, b);
}

expression sl::max(field a, field b)
{
    expression expra = a;
    expression exprb = b;
    return max(expra, exprb);
}

expression sl::max(parameter a, parameter b)
{
    expression expra = a;
    expression exprb = b;
    return max(expra, exprb);
}

expression sl::min(expression a, expression b)
{
    a.reuseit(); b.reuseit();
    
    return ifpositive(a-b, b, a);
}

expression sl::min(field a, field b)
{
    expression expra = a;
    expression exprb = b;
    return min(expra, exprb);
}

expression sl::min(parameter a, parameter b)
{
    expression expra = a;
    expression exprb = b;
    return min(expra, exprb);
}


expression sl::on(int physreg, expression expr, bool errorifnotfound)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return expr.on(physreg, NULL, errorifnotfound);
}

expression sl::on(int physreg, expression coordshift, expression expr, bool errorifnotfound)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return expr.on(physreg, &coordshift, errorifnotfound);
}

expression sl::comp(int selectedcomp, expression input)
{
    std::vector<expression> mycomp(input.countcolumns());
    for (int i = 0; i < input.countcolumns(); i++)
        mycomp[i] = input.at(selectedcomp,i);
    return expression(1, input.countcolumns(), mycomp);
}

expression sl::compx(expression input) { return comp(0,input); }
expression sl::compy(expression input) { return comp(1,input); }
expression sl::compz(expression input) { return comp(2,input); }

expression sl::entry(int row, int col, expression input) { return input.at(row,col); }

expression sl::eye(int size)
{
    if (size < 0)
    {
        std::cout << "Error in 'sl' namespace: cannot create a " << size << "x" << size << " identity matrix" << std::endl;
        abort();
    }

    std::vector<expression> exprs(size);
    for (int i = 0; i < size; i++)
        exprs[i] = 1.0;
        
    expression output(size, size, exprs);
    
    return output;
}

expression sl::transpose(expression input) { return input.transpose(); }
expression sl::inverse(expression input) { return input.invert(); }
expression sl::determinant(expression input) { return input.determinant(); }

expression sl::grad(expression input)
{
    if (input.countcolumns() != 1 || input.countrows() > 3)
    {
        std::cout << "Error in 'sl' namespace: can only take the gradient of a scalar or an up to length 3 column vector" << std::endl;
        abort();
    }

    // Cylindrical transformation of the gradient of a vector (different than of a scalar):
    if (universe::isaxisymmetric && input.countrows() > 1)
    {
        // Coordinate x on the deformed mesh:
        expression x = input.jac(2,2);
        // The output must be a nx3 matrix:
        if (input.countrows() == 2)
            return expression(2,3, {compx(dx(input)),compx(dy(input)),0, compy(dx(input)),compy(dy(input)),0});
        if (input.countrows() == 3)
            return expression(3,3, {compx(dx(input)),compx(dy(input)),-1.0/x*compz(input), compy(dx(input)),compy(dy(input)),0, compz(dx(input)),compz(dy(input)),1.0/x*compx(input)});
    }

    int problemdimension = universe::getrawmesh()->getmeshdimension();
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

expression sl::div(expression input)
{
    if (input.countcolumns() != 1 || input.countrows() > 3)
    {
        std::cout << "Error in 'sl' namespace: can only take the divergence of an up to length 3 column vector" << std::endl;
        abort();
    }

    if (universe::isaxisymmetric)
    {
        // Coordinate x on the deformed mesh:
        expression x = input.jac(2,2);
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

    abort(); // fix return warning
}

expression sl::curl(expression input)
{
    bool ishcurlfield = false;
    std::vector<std::pair<std::string,expression>> inrefcoord = input.getinrefcoord();
    if (inrefcoord.size() > 0)
    {
        ishcurlfield = ((inrefcoord[0].first) == "hcurl");
        input = inrefcoord[0].second;
    }

    if (input.countcolumns() > 1 || input.countrows() > 3)
    {
        std::cout << "Error in 'sl' namespace: can only take the curl of an up to length 3 column vector" << std::endl;
        abort();
    }

    // The curl of a hcurl type field is computed in a special way:
    if (ishcurlfield == false)
    {
        if (universe::isaxisymmetric)
        {
            // Coordinate x on the deformed mesh:
            expression x = input.jac(2,2);
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
    
    abort(); // fix return warning
}

expression sl::crossproduct(expression a, expression b)
{
    if (a.countcolumns() != 1 || b.countcolumns() != 1 || a.countrows() > 3 || b.countrows() > 3)
    {
        std::cout << "Error in 'sl' namespace: can only take the cross product of up to length 3 column vectors" << std::endl;
        abort();
    }

    a = a.resize(3,1);
    b = b.resize(3,1);

    expression a1 = compx(a), a2 = compy(a), a3 = compz(a);
    expression b1 = compx(b), b2 = compy(b), b3 = compz(b);

    expression crossprodexpr = expression(3,1, { a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1 });

    return crossprodexpr;
}

expression sl::doubledotproduct(expression a, expression b)
{
    if (a.countcolumns() != b.countcolumns() || a.countrows() != b.countrows())
    {
        std::cout << "Error in 'sl' namespace: dimension mismatch for double dot product" << std::endl;
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

expression sl::trace(expression a)
{
    if (a.countcolumns() != a.countrows())
    {
        std::cout << "Error in 'sl' namespace: can only get the trace of a square matrix" << std::endl;
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

std::vector<expression> sl::rotation(double alphax, double alphay, double alphaz, std::string type)
{
    double pi = getpi();
    double ax = alphax*pi/180.0;
    double ay = alphay*pi/180.0;
    double az = alphaz*pi/180.0;
    
    double tx,ty,tz,c,s;

    if (type == "")
    {    
        tx = ax; ty = ay; tz = az;
    
        c = std::cos(tx); s = std::sin(tx);
        densemat Rx(3,3, { 1,0,0, 0,c,-s, 0,s,c });
        c = std::cos(ty); s = std::sin(ty);
        densemat Ry(3,3, { c,0,s, 0,1,0, -s,0,c });
        c = std::cos(tz); s = std::sin(tz);
        densemat Rz(3,3, { c,-s,0, s,c,0, 0,0,1 });
        
        densemat R = Rz.multiply(Ry.multiply(Rx));
        
        double* Rval = R.getvalues();
        
        std::vector<expression> exprs(9);
        for (int i = 0; i < 9; i++)
        {
            if (std::abs(Rval[i]) < 1e-12)
                Rval[i] = 0;
            exprs[i] = expression(Rval[i]);
        }
        expression Rexpr(3,3, exprs);
        
        return std::vector<expression>{Rexpr, transpose(Rexpr)};
    }
    if (type == "voigt")
    {
        tx = ax; ty = ay; tz = az;
    
        c = std::cos(tx); s = std::sin(tx);
        densemat Kx(6,6, { 1,0,0,0,0,0, 0,c*c,s*s,-2.0*c*s,0,0, 0,s*s,c*c,2.0*c*s,0,0, 0,c*s,-c*s,c*c-s*s,0,0, 0,0,0,0,c,s, 0,0,0,0,-s,c });
        c = std::cos(ty); s = std::sin(ty);
        densemat Ky(6,6, { c*c,0,s*s,0,2.0*c*s,0, 0,1,0,0,0,0, s*s,0,c*c,0,-2.0*c*s,0, 0,0,0,c,0,-s, -c*s,0,c*s,0,c*c-s*s,0, 0,0,0,s,0,c });
        c = std::cos(tz); s = std::sin(tz);
        densemat Kz(6,6, { c*c,s*s,0,0,0,-2.0*c*s, s*s,c*c,0,0,0,2.0*c*s, 0,0,1,0,0,0, 0,0,0,c,s,0, 0,0,0,-s,c,0, c*s,-c*s,0,0,0,c*c-s*s });

        densemat K = Kz.multiply(Ky.multiply(Kx));
        
        double* Kval = K.getvalues();
        
        std::vector<expression> exprs(36);
        for (int i = 0; i < 36; i++)
        {
            if (std::abs(Kval[i]) < 1e-12)
                Kval[i] = 0;
            exprs[i] = expression(Kval[i]);
        }
        expression Kexpr(6,6, exprs);
        
        
        tx = -ax; ty = -ay; tz = -az;
    
        c = std::cos(tx); s = std::sin(tx);
        densemat invKx(6,6, { 1,0,0,0,0,0, 0,c*c,s*s,-2.0*c*s,0,0, 0,s*s,c*c,2.0*c*s,0,0, 0,c*s,-c*s,c*c-s*s,0,0, 0,0,0,0,c,s, 0,0,0,0,-s,c });
        c = std::cos(ty); s = std::sin(ty);
        densemat invKy(6,6, { c*c,0,s*s,0,2.0*c*s,0, 0,1,0,0,0,0, s*s,0,c*c,0,-2.0*c*s,0, 0,0,0,c,0,-s, -c*s,0,c*s,0,c*c-s*s,0, 0,0,0,s,0,c });
        c = std::cos(tz); s = std::sin(tz);
        densemat invKz(6,6, { c*c,s*s,0,0,0,-2.0*c*s, s*s,c*c,0,0,0,2.0*c*s, 0,0,1,0,0,0, 0,0,0,c,s,0, 0,0,0,-s,c,0, c*s,-c*s,0,0,0,c*c-s*s });

        densemat invK = invKx.multiply(invKy.multiply(invKz));
        
        double* invKval = invK.getvalues();
        
        std::vector<expression> invexprs(36);
        for (int i = 0; i < 36; i++)
        {
            if (std::abs(invKval[i]) < 1e-12)
                invKval[i] = 0;
            invexprs[i] = expression(invKval[i]);
        }
        expression invKexpr(6,6, invexprs);
        
        
        return std::vector<expression>{Kexpr, invKexpr};
    }
    
    std::cout << "Error in 'sl' namespace: rotation expected a type '' or 'voigt'" << std::endl;
    abort();
}

integration sl::integral(int physreg, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return integration(physreg, tointegrate, integrationorderdelta, blocknumber);
}

integration sl::integral(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return integration(physreg, meshdeform, tointegrate, integrationorderdelta, blocknumber);
}

integration sl::integral(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return integration(physreg, numcoefharms, tointegrate, integrationorderdelta, blocknumber);
}

integration sl::integral(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta, int blocknumber)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return integration(physreg, numcoefharms, meshdeform, tointegrate, integrationorderdelta, blocknumber);
}

expression sl::dof(expression input)
{
    return input.dof(-1);
}

expression sl::dof(expression input, int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return input.dof(physreg);
}

expression sl::tf(expression input)
{
    return input.tf(-1);
}

expression sl::tf(expression input, int physreg)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({physreg});
    return input.tf(physreg);
}

expression sl::athp(expression expr, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt)
{
    int m = expr.countrows();
    int n = expr.countcolumns();
    
    std::vector<expression> subexprs(m*n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::shared_ptr<opathp> op(new opathp(expr.getoperationinarray(i,j), rm, pt));
            subexprs[i*n+j] = expression(op);
        }
    }
    
    return expression(m, n, subexprs);
}

bool sl::adapt(int verbosity)
{
    if (universe::getrawmesh()->getdtracker()->isdefined() && slmpi::count() > 1)
    {
        std::cout << "Error in 'sl' namespace: call 'alladapt' instead of 'adapt' for multi-rank DDM" << std::endl;
        abort();
    }
    
    return universe::getrawmesh()->adapthp(verbosity);
}

bool sl::alladapt(int verbosity)
{
    return universe::getrawmesh()->adapthp(verbosity);
}

expression sl::zienkiewiczzhu(expression input)
{
    std::vector<int> alldisjregs(universe::getrawmesh()->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(input.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'sl' namespace: in 'zienkiewiczzhu' cannot have a multiharmonic expression as argument" << std::endl;
        abort();
    }

    int m = input.countrows();
    int n = input.countcolumns();

    std::vector<expression> zzexprs(m*n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (input.getoperationinarray(i,j)->isconstant())
                zzexprs[i*n+j] = 0;
            else
            {
                std::shared_ptr<opestimator> op(new opestimator("zienkiewiczzhu", input.getoperationinarray(i,j)));
                zzexprs[i*n+j] = expression(op);
            }
        }
    }
    
    return norm(expression(m, n, zzexprs));
}

expression sl::array1x1(expression term11)
{
    std::vector<expression> terms = {term11};
    return expression(1,1, terms);
}

expression sl::array1x2(expression term11, expression term12)
{
    return expression(1,2, {term11, term12});
}

expression sl::array1x3(expression term11, expression term12, expression term13)
{
    return expression(1,3, {term11, term12, term13});
}

expression sl::array2x1(expression term11, expression term21)
{
    return expression(2,1, {term11, term21});
}

expression sl::array2x2(expression term11, expression term12, expression term21, expression term22)
{
    return expression(2,2, {term11, term12, term21, term22});
}

expression sl::array2x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23)
{
    return expression(2,3, {term11, term12, term13, term21, term22, term23});
}

expression sl::array3x1(expression term11, expression term21, expression term31)
{
    return expression(3,1, {term11, term21, term31});
}

expression sl::array3x2(expression term11, expression term12, expression term21, expression term22, expression term31, expression term32)
{
    return expression(3,2, {term11, term12, term21, term22, term31, term32});
}

expression sl::array3x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23, expression term31, expression term32, expression term33)
{
    return expression(3,3, {term11, term12, term13, term21, term22, term23, term31, term32, term33});
}

vec sl::solve(mat A, vec b, std::string soltype, bool diagscaling)
{
    if (soltype != "lu" && soltype != "cholesky")
    {
        std::cout << "Error in 'sl' namespace: unknown direct solver type '" << soltype << "' (use 'lu' or 'cholesky')" << std::endl;
        abort();
    }
    if (A.countrows() != b.size())
    {
        std::cout << "Error in 'sl' namespace: direct solve of Ax = b failed (size of A and b do not match)" << std::endl;
        abort();
    }

    if (A.getpointer() == NULL || b.getpointer() == NULL)
    {
        std::cout << "Error in 'sl' namespace: direct solve of Ax = b failed (A or b is undefined)" << std::endl;
        abort();
    }
    
    vec breduced = A.eliminate(b);
    
    Vec bpetsc = breduced.getpetsc();
    Mat Apetsc = A.getapetsc();

    vec sol(std::shared_ptr<rawvec>(new rawvec(breduced.getpointer()->getdofmanager())));
    Vec solpetsc = sol.getpetsc();

    KSP* ksp = A.getpointer()->getksp();

    if (A.getpointer()->isfactored() == false)
    {
        PC pc;
        KSPCreate(PETSC_COMM_SELF, ksp);
        KSPSetOperators(*ksp, Apetsc, Apetsc);
        // Perform a diagonal scaling for improved matrix conditionning.
        // This modifies the matrix A and right handside b!
        if (diagscaling == true)
            KSPSetDiagonalScale(*ksp, PETSC_TRUE);
        KSPSetFromOptions(*ksp);

        KSPGetPC(*ksp,&pc);
        if (soltype == "lu")
            PCSetType(pc,PCLU);
        if (soltype == "cholesky")
            PCSetType(pc,PCCHOLESKY);
        PCFactorSetMatSolverType(pc, universe::solvertype);
    }

    KSPSolve(*ksp, bpetsc, solpetsc);

    A.getpointer()->isfactored(true);

    if (A.getpointer()->isfactorizationreuseallowed() == false)
    {
        KSPDestroy(ksp);
        A.getpointer()->isfactored(false);
    }

    return A.xbmerge(sol, b);
}

std::vector<vec> sl::solve(mat A, std::vector<vec> b, std::string soltype)
{
    if (soltype != "lu" && soltype != "cholesky")
    {
        std::cout << "Error in 'sl' namespace: unknown direct solver type '" << soltype << "' (use 'lu' or 'cholesky')" << std::endl;
        abort();
    }
    for (int i = 0; i < b.size(); i++)
    {
        if (A.countrows() != b[i].size())
        {
            std::cout << "Error in 'sl' namespace: multi-rhs direct solve of Ax = b failed (size of A and at least one rhs do not match)" << std::endl;
            abort();
        }
        if (A.getpointer() == NULL || b[i].getpointer() == NULL)
        {
            std::cout << "Error in 'sl' namespace: multi-rhs direct solve of Ax = b failed (A or at least one rhs is undefined)" << std::endl;
            abort();
        }
    }
    
    if (b.size() == 0)
        return {};
        
    int numrhs = b.size();
    
    std::vector<vec> breduced(numrhs);
    for (int i = 0; i < numrhs; i++)
        breduced[i] = A.eliminate(b[i]);

    int len = breduced[0].size();
    
    // Concatenate rhs vecs to densemat:
    densemat rhs(numrhs, len);
    double* rhsptr = rhs.getvalues();
    for (int i = 0; i < numrhs; i++)
    {
        densemat vecvals = breduced[i].getallvalues();
        double* vecvalsptr = vecvals.getvalues();
        for (int j = 0; j < len; j++)
            rhsptr[i*len+j] = vecvalsptr[j];
    }
    
    // Solve multi-rhs:
    densemat sols = solve(A, rhs, soltype);
    double* solsptr = sols.getvalues();

    // Extract 'sols' rows to sol vecs:
    densemat vals(len,1);
    double* valsptr = vals.getvalues();
    std::vector<vec> outvecs(numrhs);
    for (int i = 0; i < numrhs; i++)
    {
        for (int j = 0; j < len; j++)
            valsptr[j] = solsptr[i*len+j];
    
        outvecs[i] = vec(std::shared_ptr<rawvec>(new rawvec(b[i].getpointer()->getdofmanager())));
        outvecs[i].setvalues(A.getainds(), vals);
        outvecs[i].setvalues(A.getdinds(), b[i].getvalues(A.getdinds()));
    }
    
    return outvecs;
}

densemat sl::solve(mat A, densemat b, std::string soltype)
{
    int numrhs = b.countrows();
    int len = b.countcolumns();
 
    Mat Apetsc = A.getapetsc();
    
    KSP* ksp = A.getpointer()->getksp();
    PC pc;
    if (A.getpointer()->isfactored() == false)
    {
        KSPCreate(PETSC_COMM_SELF, ksp);
        KSPSetOperators(*ksp, Apetsc, Apetsc);
        KSPSetFromOptions(*ksp);

        KSPGetPC(*ksp,&pc);
        if (soltype == "lu")
            PCSetType(pc,PCLU);
        if (soltype == "cholesky")
            PCSetType(pc,PCCHOLESKY);
        PCFactorSetMatSolverType(pc, universe::solvertype);
        PCSetUp(pc);
    }
    else
        KSPGetPC(*ksp,&pc);
        
    PCFactorGetMatrix(pc, &Apetsc);

    densemat densesols(numrhs, len);

    Mat sols, rhses;
    MatCreateSeqDense(PETSC_COMM_SELF, len, numrhs, densesols.getvalues(), &sols);
    MatCreateSeqDense(PETSC_COMM_SELF, len, numrhs, b.getvalues(), &rhses);
    
    // 'rhs' and 'sols' are considered column major in petsc.
    MatMatSolve(Apetsc, rhses, sols);
    
    MatDestroy(&sols);
    MatDestroy(&rhses);

    A.getpointer()->isfactored(true);

    if (A.getpointer()->isfactorizationreuseallowed() == false)
    {
        KSPDestroy(ksp);
        A.getpointer()->isfactored(false);
    }
    
    return densesols;
}

int mykspmonitor(KSP ksp, PetscInt iter, PetscReal resnorm, void* unused)
{
    std::cout << iter << " KSP residual norm " << resnorm << std::endl;
    return 0;
}

void sl::solve(mat A, vec b, vec sol, double& relrestol, int& maxnumit, std::string soltype, std::string precondtype, int verbosity, bool diagscaling)
{
    if (soltype != "gmres" && soltype != "bicgstab")
    {
        std::cout << "Error in 'sl' namespace: unknown iterative solver type '" << soltype << "' (use 'gmres' or 'bicgstab')" << std::endl;
        abort();
    }
    if (precondtype != "ilu" && precondtype != "sor" && precondtype != "gamg")
    {
        std::cout << "Error in 'sl' namespace: unknown preconditioner type '" << precondtype << "' (use 'ilu', 'sor' or 'gamg')" << std::endl;
        abort();
    }
    if (A.countrows() != b.size())
    {
        std::cout << "Error in 'sl' namespace: iterative solve of Ax = b failed (size of A and b do not match)" << std::endl;
        abort();
    }
    if (A.countrows() != sol.size())
    {
        std::cout << "Error in 'sl' namespace: iterative solve of Ax = b failed (size of A and x do not match)" << std::endl;
        abort();
    }

    if (A.getpointer() == NULL || b.getpointer() == NULL || sol.getpointer() == NULL)
    {
        std::cout << "Error in 'sl' namespace: iterative solve of Ax = b failed (A, x or b is undefined)" << std::endl;
        abort();
    }
    
    vec breduced = A.eliminate(b);
    vec sola = sol.extract(A.getainds());

    Vec bpetsc = breduced.getpetsc();
    Mat Apetsc = A.getapetsc();

    Vec solpetsc = sola.getpetsc();

    KSP* ksp = A.getpointer()->getksp();

    KSPCreate(PETSC_COMM_SELF, ksp);
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
    if (precondtype == "gamg")
        PCSetType(pc,PCGAMG);

    KSPSolve(*ksp, bpetsc, solpetsc);

    // Get the number of required iterations and the residual norm:
    KSPGetIterationNumber(*ksp, &maxnumit);
    KSPGetResidualNorm(*ksp, &relrestol);

    KSPDestroy(ksp);
    
    sol.setvalues(A.getainds(), sola.getallvalues());
    sol.setvalues(A.getdinds(), b.getvalues(A.getdinds()));
}

void sl::exchange(std::vector<int> targetranks, std::vector<densemat> sends, std::vector<densemat> receives)
{
    int numtargets = targetranks.size();

    if (numtargets == 0)
        return;

    std::vector<int> sendlens(numtargets), reclens(numtargets);
    std::vector<double*> sendbuffers(numtargets), recbuffers(numtargets);

    for (int i = 0; i < numtargets; i++)
    {
        sendlens[i] = sends[i].count();
        reclens[i] = receives[i].count();

        sendbuffers[i] = sends[i].getvalues();
        recbuffers[i] = receives[i].getvalues();
    }

    slmpi::exchange(targetranks, sendlens, sendbuffers, reclens, recbuffers);
}

std::vector<double> sl::gmres(densemat (*mymatmult)(densemat), densemat b, densemat x, double relrestol, int maxnumit, int verbosity)
{   
    if (b.countrows() != x.countrows() || b.countcolumns() != 1 || x.countcolumns() != 1)
    {
        std::cout << "Error in 'sl' namespace: in function gmres expected a column vector of same size for b and x" << std::endl;
        abort();
    }

    // Fragment size:
    int n = b.count();
    
    // Initialize the 1D vectors:
    std::vector<double> sn(maxnumit, 0.0);
    std::vector<double> cs(maxnumit, 0.0);
    std::vector<double> beta(maxnumit+1, 0.0);
    std::vector<double> relresvec(maxnumit+1, 0.0);
    
    double* xptr = x.getvalues();
    double* bptr = b.getvalues();

    // Compute r = b - A * x and the initial relative residual:
    double normb = 0.0; double normr = 0.0;
    densemat r = mymatmult(x.copy());
    double* rptr = r.getvalues();
    
    if (r.countrows() != n || r.countcolumns() != 1)
    {
        std::cout << "Error in 'sl' namespace: in function gmres the matrix product function call returned a densemat of wrong size on rank " << slmpi::getrank() << std::endl;
        abort();
    }
    
    for (int i = 0; i < n; i++)
    {
        rptr[i] = bptr[i] - rptr[i];
        normb += bptr[i]*bptr[i];
        normr += rptr[i]*rptr[i];
    }
    
    std::vector<double> norms = {normb, normr};
    
    slmpi::sum(norms); // reduce on all ranks
    
    normb = std::sqrt(norms[0]);
    normr = std::sqrt(norms[1]);
    beta[0] = normr;
    
    // All zero solution in case b is all zero:
    if (normb == 0)
    {
        for (int i = 0; i < n; i++)
            xptr[i] = 0.0;
        return {0.0};
    }
    
    relresvec[0] = normr/normb;
        
    // Holder for all Krylov vectors (one on each row):
    densemat Q(maxnumit+1, n);
    double* Qptr = Q.getvalues();
    
    // First row is the normed residual:
    double invnormr = 1.0/normr;
    for (int i = 0; i < n; i++)
        Qptr[i] = invnormr * rptr[i];
        
    // Hessenberg matrix (columnwise upper triangular {r0c0,r0c1,r1c1,r0c2,...}):
    densemat H(1, ((1+maxnumit)*maxnumit)/2 + 1, 0.0); // +1 because arnoldi returns length k+2
    double* Hptr = H.getvalues();
    
    // GMRES iteration:
    int k = 0;
    for (k = 0; k < maxnumit; k++)
    {
        if (verbosity > 0)
            std::cout << "gmres @" << k << " -> " << relresvec[k] << std::endl;
            
        if (relresvec[k] <= relrestol)
            break;
            
        // Run Arnoldi:
        std::vector<double> h = gentools::arnoldi(mymatmult, Q, k);
        
        // Write h to H (can exceed by 1 the reduced Hessenberg matrix size):
        for (int i = 0; i < h.size(); i++)
            Hptr[((1+k)*k)/2 + i] = h[i];

        // Eliminate the last element in the kth column of H and update the rotation matrix:
        gentools::applygivensrotation(Hptr+((1+k)*k)/2, cs, sn, k);
        
        // Update the residual vector:
        beta[k+1] = -sn[k] * beta[k];
        beta[k] = cs[k] * beta[k];
        
        relresvec[k+1] = std::abs(beta[k+1]) / normb;
    }

    if (k > 0)
    {
        // Calculate the solution:
        densemat y(1,k);
        gentools::solveuppertriangular(k, Hptr, &beta[0], y.getvalues());
        densemat Qy = y.multiply(Q.getresized(k,n));
        x.add(Qy);
    }

    relresvec.resize(k+1);
    
    return relresvec;
}

void sl::mapdofs(std::shared_ptr<dofmanager> dm, std::vector<std::shared_ptr<rawfield>> rfs, std::vector<bool> isdimactive, std::vector<indexmat>& sendinds, std::vector<indexmat>& recvinds)
{
    std::shared_ptr<dtracker> dt = universe::getrawmesh()->getdtracker();
    
    elements* els = dt->getrawmesh()->getelements();
    disjointregions* drs = dt->getrawmesh()->getdisjointregions();
    physicalregions* prs = dt->getrawmesh()->getphysicalregions();
    
    int rank = slmpi::getrank();
    
    std::vector<std::vector<std::vector<int>>>* elemmap = dt->getmap();
    
    int numrawfields = rfs.size();
    int numdisjregs = drs->count();
    int numneighbours = dt->countneighbours();
    
    // Get the disjoint regions in each overlap/no-overlap interface:
    std::vector<std::vector<bool>> isdisjregininnerinterface(numneighbours, std::vector<bool>(numdisjregs, false));
    std::vector<std::vector<bool>> isdisjreginouterinterface(numneighbours, std::vector<bool>(numdisjregs, false));
    for (int n = 0; n < numneighbours; n++)
    {
        int cn = dt->getneighbour(n);
        
        if (dt->isoverlap())
        {
            for (int dim = 0; dim < 3; dim++)
            {
                if (isdimactive[dim])
                {
                    int cr = dt->getinneroverlapinterface(cn, dim);
                    if (cr >= 0)
                    {
                        std::vector<int> cdrs = prs->get(cr)->getdisjointregions(-1);
                        for (int i = 0; i < cdrs.size(); i++)
                            isdisjregininnerinterface[n][cdrs[i]] = true;
                    }
                    cr = dt->getouteroverlapinterface(cn, dim);
                    if (cr >= 0)
                    {
                        std::vector<int> cdrs = prs->get(cr)->getdisjointregions(-1);
                        for (int i = 0; i < cdrs.size(); i++)
                            isdisjreginouterinterface[n][cdrs[i]] = true;
                    }
                }
            }
        }
        else
        {
            for (int dim = 0; dim < 3; dim++)
            {
                if (isdimactive[dim])
                {
                    int cr = dt->getnooverlapinterface(cn, dim);
                    if (cr >= 0)
                    {
                        std::vector<int> cdrs = prs->get(cr)->getdisjointregions(-1);
                        for (int i = 0; i < cdrs.size(); i++)
                            isdisjregininnerinterface[n][cdrs[i]] = true;
                    }
                }
            }
            isdisjreginouterinterface[n] = isdisjregininnerinterface[n];
        }
    }
    
    // Count the number of dofs to send for each rawfield and the number of data send ranges:
    std::vector<std::vector<int>> numsenddofsperfield(numneighbours, std::vector<int>(numrawfields, 0));
    std::vector<std::vector<int>> numsendrangesperfield(numneighbours, std::vector<int>(numrawfields, 0));
    std::vector<std::vector<int>> numexpectedrecvdofsperfield(numneighbours, std::vector<int>(numrawfields, 0));
    for (int n = 0; n < numneighbours; n++)
    {
        for (int r = 0; r < numrawfields; r++)
        {
            dm->selectfield(rfs[r]);
            for (int d = 0; d < numdisjregs; d++)
            {
                int ne = drs->countelements(d);
                int nff = dm->countformfunctions(d);

                if (isdisjregininnerinterface[n][d])
                {
                    numsenddofsperfield[n][r] += ne*nff;
                    if (nff > 0)
                        numsendrangesperfield[n][r]++;
                }
                if (isdisjreginouterinterface[n][d])
                    numexpectedrecvdofsperfield[n][r] += ne*nff;
            }
        }
    }
    
    // Preallocate the outputs:
    sendinds = std::vector<indexmat>(numneighbours);
    recvinds = std::vector<indexmat>(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        sendinds[n] = indexmat(gentools::sum(numsenddofsperfield[n]), 1);
        recvinds[n] = indexmat(gentools::sum(numexpectedrecvdofsperfield[n]), 1);
    }
    
    // Create the send range data and populate the send indexes:
    std::vector<std::vector<int>> dofrangesforneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
    {
        int index = 0;
        int* sivals = sendinds[n].getvalues();
        int totnumsendranges = gentools::sum(numsendrangesperfield[n]);
        
        // Format is {numrfs, numsenddofsrf0,...,numsenddofsrfn, numrangesrf0,...,numrangesrfn, type0,rbe0,ree0,nff0, type1,...}:    
        dofrangesforneighbours[n] = std::vector<int>(1 + 2*numrawfields + 4*totnumsendranges);
        dofrangesforneighbours[n][0] = numrawfields;
        
        int ind = 1+2*numrawfields;
        for (int r = 0; r < numrawfields; r++)
        {
            dofrangesforneighbours[n][1+r] = numsenddofsperfield[n][r];
            dofrangesforneighbours[n][1+numrawfields+r] = numsendrangesperfield[n][r];
            
            dm->selectfield(rfs[r]);
            for (int d = 0; d < numdisjregs; d++)
            {
                if (isdisjregininnerinterface[n][d])
                {
                    int elemtype = drs->getelementtypenumber(d);
                    int rbe = drs->getrangebegin(d);
                    int ree = drs->getrangeend(d);
                    int ne = ree-rbe+1;
                    int nff = dm->countformfunctions(d);
                    if (nff > 0)
                    {
                        dofrangesforneighbours[n][ind+0] = elemtype;
                        dofrangesforneighbours[n][ind+1] = rbe;
                        dofrangesforneighbours[n][ind+2] = ree;
                        dofrangesforneighbours[n][ind+3] = nff;
                        
                        for (int f = 0; f < nff; f++)
                        {
                            int rb = dm->getrangebegin(d, f);
                            
                            for (int e = 0; e < ne; e++)
                                sivals[index+e] = rb+e;
                            index += ne;
                        }
                        ind += 4;
                    }
                }
            }
        }
    }
    
    std::vector<int> sendlens(numneighbours), recvlens;
    for (int n = 0; n < numneighbours; n++)
        sendlens[n] = dofrangesforneighbours[n].size();
        
    slmpi::exchange(dt->getneighbours(), sendlens, recvlens);
    
    std::vector<std::vector<int>> dofrangesfromneighbours(numneighbours);
    for (int n = 0; n < numneighbours; n++)
        dofrangesfromneighbours[n].resize(recvlens[n]);
        
    slmpi::exchange(dt->getneighbours(), dofrangesforneighbours, dofrangesfromneighbours);
    
    // Sanity check. The number of rawfields is allowed to change from a domain to another.
    for (int n = 0; n < numneighbours; n++)
    {
        int numrecvrawfields = dofrangesfromneighbours[n][0];
        
        std::vector<int> numdofsexpected = {};
        for (int r = 0; r < numrawfields; r++)
        {
            if (numexpectedrecvdofsperfield[n][r] > 0)
                numdofsexpected.push_back(numexpectedrecvdofsperfield[n][r]);
        }
        std::vector<int> numdofsrecv = {};
        for (int r = 0; r < numrecvrawfields; r++)
        {
            if (dofrangesfromneighbours[n][1+r] > 0)
                numdofsrecv.push_back(dofrangesfromneighbours[n][1+r]);
        }
        if (numdofsexpected != numdofsrecv)
        {
            std::cout << "Error in 'sl' namespace: DDM interface data size does not match for at least one field between ranks " << rank << " and " << dt->getneighbour(n) << std::endl;
            std::cout << "Make sure that the fields are provided in the same order to the dof manager and that the interpolation orders match on the DDM interfaces" << std::endl; 
            abort();
        }
    }

    // Populate the receive indexes (already preallocated above):
    for (int n = 0; n < numneighbours; n++)
    {
        int numrecvrawfields = dofrangesfromneighbours[n][0];
        
        int index = 0;
        int* rivals = recvinds[n].getvalues();
        
        int ind = 1+2*numrecvrawfields;
        
        int recvrfindex = 0;
        for (int r = 0; r < numrawfields; r++)
        {
            if (numexpectedrecvdofsperfield[n][r] == 0)
                continue;

            while (dofrangesfromneighbours[n][1+recvrfindex] == 0)
                recvrfindex++;
            
            dm->selectfield(rfs[r]);
                
            int recvnumranges = dofrangesfromneighbours[n][1+numrecvrawfields+recvrfindex];
            for (int i = 0; i < recvnumranges; i++)
            {
                int elemtype = dofrangesfromneighbours[n][ind+0];
                int rbe = dofrangesfromneighbours[n][ind+1];
                int ree = dofrangesfromneighbours[n][ind+2];
                int nff = dofrangesfromneighbours[n][ind+3];
                
                int ne = ree-rbe+1;
                
                for (int f = 0; f < nff; f++)
                {
                    for (int e = 0; e < ne; e++)
                    {
                        int ce = elemmap->at(n)[elemtype][rbe+e];
                        int cdr = els->getdisjointregion(elemtype, ce);
                        int crbe = drs->getrangebegin(cdr);
                        int crb = dm->getrangebegin(cdr, f);
                        
                        rivals[index+e] = crb + ce-crbe;
                    }
                    index += ne;
                }
                ind += 4;
            }
            recvrfindex++;            
        }
    }
}

std::vector<double> sl::linspace(double a, double b, int num)
{
    if (num < 0)
    {
        std::cout << "Error in 'sl' namespace: cannot call 'linspace' for " << num << " points" << std::endl;
        abort();
    }
    if (num == 0)
        return {};
    if (num == 1)
        return {b};

    double step = (b-a)/(num-1.0);

    std::vector<double> output(num);
    for (int i = 0; i < num; i++)
        output[i] = a + step*i;    

    return output;
}

std::vector<double> sl::logspace(double a, double b, int num, double basis)
{
    if (num < 0)
    {
        std::cout << "Error in 'sl' namespace: cannot call 'logspace' for " << num << " points" << std::endl;
        abort();
    }
    if (num == 0)
        return {};
    if (num == 1)
        return {std::pow(basis, b)};

    double step = (b-a)/(num-1.0);

    std::vector<double> output(num);
    for (int i = 0; i < num; i++)
        output[i] = std::pow(basis, a + step*i);    

    return output;
}

expression sl::dbtoneper(expression toconvert)
{
    if (toconvert.iszero())
        return toconvert;

    return ( 0.1151292546497023 * toconvert );
}

void sl::setdata(vec invec)
{
    if (invec.getpointer() == NULL)
        return;
    
    // Get all fields in the vec structure:
    std::vector<std::shared_ptr<rawfield>> allfields = invec.getpointer()->getdofmanager()->getfields();
    
    for (int i = 0; i < allfields.size(); i++)
        allfields[i]->setdata(-1, invec|field(allfields[i]));
        
    invec.getpointer()->setvaluestoports();
}


////////// PREDEFINED OPERATORS

expression sl::strain(expression input)
{
    if ((input.countrows() != 2 && input.countrows() != 3) || (input.countcolumns() != 1 && input.countrows() != input.countcolumns()))
    {
        std::cout << "Error in 'sl' namespace: can only compute the strains of a 2x1 or 3x1 column vector or its gradient" << std::endl;
        abort();
    }

    expression gradu = input;
    if (input.countcolumns() == 1)
        gradu = sl::grad(input);

    if (input.countrows() == 2)
        return expression(3,1,{gradu.at(0,0), gradu.at(1,1), gradu.at(1,0) + gradu.at(0,1)});
    if (input.countrows() == 3)
        return expression(6,1,{gradu.at(0,0), gradu.at(1,1), gradu.at(2,2), gradu.at(2,1) + gradu.at(1,2), gradu.at(0,2) + gradu.at(2,0), gradu.at(0,1) + gradu.at(1,0)});
        
    abort(); // fix return warning
}

expression sl::greenlagrangestrain(expression input)
{
    if ((input.countrows() != 2 && input.countrows() != 3) || (input.countcolumns() != 1 && input.countrows() != input.countcolumns()))
    {
        std::cout << "Error in 'sl' namespace: can only compute the green-lagrange strains of a 2x1 or 3x1 column vector or its gradient" << std::endl;
        abort();
    }

    expression gradu = input;
    if (input.countcolumns() == 1)
        gradu = sl::grad(input);

    // This can be called since gradu is nonlinear in u and can thus not include a dof or tf:
    gradu.reuseit();

    if (input.countrows() == 2)
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
    if (input.countrows() == 3)
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
    
    abort(); // fix return warning
}

expression sl::vonmises(expression stress)
{
    if (stress.countcolumns() != 1 || stress.countrows() != 6)
    {
        std::cout << "Error in 'sl' namespace: expected the 3D stress tensor in Voigt notation (6 rows, 1 column)" << std::endl;
        abort();
    }

    expression s11 = stress.at(0,0), s22 = stress.at(1,0), s33 = stress.at(2,0), s23 = stress.at(3,0), s13 = stress.at(4,0), s12 = stress.at(5,0);
    s11.reuseit(); s22.reuseit(); s33.reuseit();

    return sqrt( 0.5*( pow(s11-s22,2)+pow(s22-s33,2)+pow(s33-s11,2) ) + 3.0*( pow(s12,2)+pow(s23,2)+pow(s13,2) ) );
}

expression sl::predefinedmassconservation(expression dofv, expression tfp, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant)
{
    if (isdensityconstant)
        return div(dofv)*tfp;

    if (includetimederivs)
        return ( rho*div(dofv)*tfp + dofv*gradrho*tfp + dtrho*tfp );
    else
        return ( rho*div(dofv)*tfp + dofv*gradrho*tfp );
}

expression sl::predefinedinertialforce(expression dofv, expression tfv, expression v, expression rho)
{
    rho.reuseit(); v.reuseit();

    return ( -rho*( grad(v)*dofv + grad(dofv)*v - grad(v)*v )*tfv );
}

expression sl::predefinedviscousforce(expression dofv, expression tfv, expression mu, bool isdensityconstant, bool isviscosityconstant)
{
    mu.reuseit();

    if (isdensityconstant && isviscosityconstant)
        return ( - mu*doubledotproduct(grad(dofv), grad(tfv)) );

    if (isdensityconstant)
        return ( - mu*doubledotproduct(grad(dofv), grad(tfv)) - mu*doubledotproduct(transpose(grad(dofv)), grad(tfv)) );
    else
        return ( - mu*doubledotproduct(grad(dofv), grad(tfv)) - mu*doubledotproduct(transpose(grad(dofv)), grad(tfv)) + (2.0/3.0)*mu*div(dofv)*trace(grad(tfv)) );
}


////////// PREDEFINED FORMULATIONS

std::vector<integration> sl::continuitycondition(int gamma1, int gamma2, field u1, field u2, int lagmultorder, bool errorifnotfound)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({gamma1, gamma2});
    
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    int gamma1dim = universe::getrawmesh()->getphysicalregions()->get(gamma1)->getelementdimension();
    int gamma2dim = universe::getrawmesh()->getphysicalregions()->get(gamma2)->getelementdimension();
    if (gamma1dim != gamma2dim || gamma1dim >= problemdimension)
    {
        std::cout << "Error in 'sl' namespace: expected boundary regions for gamma1 and gamma2 in 'continuitycondition'" << std::endl;
        abort();
    }
    
    std::shared_ptr<rawfield> ptr1 = u1.getpointer();
    std::shared_ptr<rawfield> ptr2 = u2.getpointer();
    
    // Make sure the fields are similar:
    if (ptr1->gettypename(false) != ptr2->gettypename(false) || ptr1->getharmonics() != ptr2->getharmonics())
    {
        std::cout << "Error in 'sl' namespace: in 'continuitycondition' expected two fields of same type and harmonic content" << std::endl;
        abort();
    }
    
    // Create the Lagrange multiplier field:
    field lambda(ptr1->gettypename(false), ptr1->getharmonics());
    lambda.noautomaticupdate();
    lambda.getpointer()->setorder(gamma1, lagmultorder, false);
    lambda.getpointer()->setorder(gamma2, lagmultorder, false);

    // Create the integration object to output:
    std::vector<integration> output(3);
    
    output[0] = integral(gamma1, dof(lambda)*tf(u1));
    output[1] = integral(gamma2, -on(gamma1, dof(lambda), errorifnotfound) * tf(u2));
    output[2] = integral(gamma1, (dof(u1) - on(gamma2, dof(u2), errorifnotfound)) * tf(lambda));

    return output;
}

std::vector<integration> sl::continuitycondition(int gamma1, int gamma2, field u1, field u2, std::vector<double> rotcent, double rotangz, double angzmod, double factor, int lagmultorder)
{       
    universe::getrawmesh()->getphysicalregions()->errorundefined({gamma1, gamma2});
    
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    int gamma1dim = universe::getrawmesh()->getphysicalregions()->get(gamma1)->getelementdimension();
    int gamma2dim = universe::getrawmesh()->getphysicalregions()->get(gamma2)->getelementdimension();
    if (gamma1dim != gamma2dim || gamma1dim >= problemdimension)
    {
        std::cout << "Error in 'sl' namespace: expected boundary regions for gamma1 and gamma2 in 'continuitycondition'" << std::endl;
        abort();
    }

    std::shared_ptr<rawfield> ptr1 = u1.getpointer();
    std::shared_ptr<rawfield> ptr2 = u2.getpointer();
    
    // Make sure the fields are similar:
    if (ptr1->gettypename(false) != ptr2->gettypename(false) || ptr1->getharmonics() != ptr2->getharmonics())
    {
        std::cout << "Error in 'sl' namespace: in 'continuitycondition' expected two fields of same type and harmonic content" << std::endl;
        abort();
    }
    
    if (rotcent.size() != 3)
    {
        std::cout << "Error in 'sl' namespace: in 'continuitycondition' expected a vector of length 3 as fifth argument" << std::endl;
        abort();
    }
    
    if (factor != -1 && factor != 1)
    {
        std::cout << "Error in 'sl' namespace: in 'continuitycondition' the factor must be -1 or 1" << std::endl;
        abort();
    }

    if (angzmod < 0.0 || angzmod > 180.0)
    {
        std::cout << "Error in 'sl' namespace: in 'continuitycondition' the angular modulo should be in range [0,180]" << std::endl;
        abort();
    }

    if (rotangz < 0.0 || rotangz > angzmod)
    {
        std::cout << "Error in 'sl' namespace: in 'continuitycondition' the rotation angle should be in range [0,angzmod]" << std::endl;
        abort();
    }
    
    // Create the Lagrange multiplier field:
    field lambda(ptr1->gettypename(false), ptr1->getharmonics());
    lambda.noautomaticupdate();
    lambda.getpointer()->setorder(gamma1, lagmultorder, false);
    lambda.getpointer()->setorder(gamma2, lagmultorder, false);
    
    int numcomp = expression(lambda).countrows();

    // Create the face to face mapping expression:
    expression mapexpr, invmapexpr;
    expression tfu = tf(u2);
    expression dofu = dof(u2);
    
    if (numcomp > 1)
    {
        tfu = tfu.resize(3,1);
        dofu = dofu.resize(3,1);
    }
    
    field x("x"), y("y"), z("z");
    
    expression centered = array3x1(x-rotcent[0],y-rotcent[1],z);
    
    expression radius = sqrt( compx(centered)*compx(centered) + compy(centered)*compy(centered) );
    
    mapexpr = centered;
    // Calculate the angle for the gamma1 coordinates:
    expression mapangle = acos( compx(mapexpr)/radius );
    mapangle = ifpositive(compy(mapexpr), mapangle, 2.0*getpi()-mapangle);
    // Take the angular modulo:
    expression mapmod = mapexpr;
    mapmod.rotate(0,0,angzmod);
    mapexpr = ifpositive(rotangz*getpi()/180.0 - mapangle, mapmod, mapexpr);
    expression doffact = ifpositive(rotangz*getpi()/180.0 - mapangle, factor, 1.0);
    mapexpr = mapexpr + array3x1(rotcent[0],rotcent[1],0) - array3x1(x,y,z);
    
    invmapexpr = centered;
    // Calculate the angle for the gamma2 coordinates:
    expression invmapangle = acos( compx(invmapexpr)/radius );
    invmapangle = ifpositive(compy(invmapexpr), invmapangle, 2.0*getpi()-invmapangle);
    // Take the angular modulo:
    expression invmapmod = invmapexpr;
    invmapmod.rotate(0,0,-angzmod);
    invmapexpr = ifpositive(invmapangle - angzmod*getpi()/180.0, invmapmod, invmapexpr);
    expression tffact = ifpositive(invmapangle - angzmod*getpi()/180.0, factor, 1.0);
    invmapexpr = invmapexpr + array3x1(rotcent[0],rotcent[1],0) - array3x1(x,y,z);
    
    if (numcomp > 1)
    {
        expression theta = ifpositive(invmapangle - angzmod*getpi()/180.0, -angzmod*getpi()/180.0, 0.0);
        expression rotmat = array3x3(cos(theta),-sin(theta),0, sin(theta),cos(theta),0, 0,0,1);
        dofu = rotmat*dofu;
        tfu = rotmat*tfu;
    }
    
    if (numcomp > 1)
    {
        tfu = tfu.resize(numcomp,1);
        dofu = dofu.resize(numcomp,1);
    }
    
    
    // Create the integration object to output:
    std::vector<integration> output(3);
    
    output[0] = integral(gamma1, dof(lambda)*tf(u1));
    output[1] = integral(gamma2, -on(gamma1, invmapexpr, dof(lambda)) * tffact * tfu);
    output[2] = integral(gamma1, (dof(u1) - doffact * on(gamma2, mapexpr, dofu)) * tf(lambda));

    return output;
}

std::vector<integration> sl::periodicitycondition(int gamma1, int gamma2, field u, std::vector<double> dat1, std::vector<double> dat2, double factor, int lagmultorder)
{
    universe::getrawmesh()->getphysicalregions()->errorundefined({gamma1, gamma2});
    
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    int gamma1dim = universe::getrawmesh()->getphysicalregions()->get(gamma1)->getelementdimension();
    int gamma2dim = universe::getrawmesh()->getphysicalregions()->get(gamma2)->getelementdimension();
    if (gamma1dim != gamma2dim || gamma1dim >= problemdimension)
    {
        std::cout << "Error in 'sl' namespace: expected boundary regions for gamma1 and gamma2 in 'periodicitycondition'" << std::endl;
        abort();
    }
    
    if (dat1.size() != 3)
    {
        std::cout << "Error in 'sl' namespace: in 'periodicitycondition' expected a vector of length 3 as fourth argument" << std::endl;
        abort();
    }
    if (dat2.size() != 1 && dat2.size() != 3)
    {
        std::cout << "Error in 'sl' namespace: in 'periodicitycondition' expected a vector of length 1 or 3 as fifth argument" << std::endl;
        abort();
    }

    // Create the Lagrange multiplier field:
    std::shared_ptr<rawfield> ptr = u.getpointer();
    field lambda(ptr->gettypename(false), ptr->getharmonics());
    lambda.noautomaticupdate();
    lambda.getpointer()->setorder(gamma1, lagmultorder, false);
    lambda.getpointer()->setorder(gamma2, lagmultorder, false);
    
    int numcomp = expression(lambda).countrows();

    // Create the face to face mapping expression:
    expression mapexpr, invmapexpr;
    expression tfu = tf(u);
    expression dofu = dof(u);
    
    if (numcomp > 1)
    {
        tfu = tfu.resize(3,1);
        dofu = dofu.resize(3,1);
    }
    
    // For a translation:
    if (dat2.size() == 1)
    {
        // First norm the direction vector:
        double normval = std::sqrt(dat1[0]*dat1[0] + dat1[1]*dat1[1] + dat1[2]*dat1[2]);
        dat1[0] = dat1[0]/normval; dat1[1] = dat1[1]/normval; dat1[2] = dat1[2]/normval;
        
        double shiftlen = dat2[0];
    
        mapexpr = array3x1(shiftlen*dat1[0], shiftlen*dat1[1], shiftlen*dat1[2]);
        invmapexpr = -mapexpr;
    }

    // For a rotation:
    if (dat2.size() == 3)
    {
        field x("x"), y("y"), z("z");
        
        mapexpr = array3x1(x-dat1[0],y-dat1[1],z-dat1[2]);
        mapexpr.rotate(dat2[0], dat2[1], dat2[2]);
        mapexpr = mapexpr + array3x1(dat1[0],dat1[1],dat1[2]) - array3x1(x,y,z);
        
        invmapexpr = array3x1(x-dat1[0],y-dat1[1],z-dat1[2]);
        invmapexpr.rotate(0,0, -dat2[2]);
        invmapexpr.rotate(0, -dat2[1], 0);
        invmapexpr.rotate(-dat2[0], 0, 0);
        invmapexpr = invmapexpr + array3x1(dat1[0],dat1[1],dat1[2]) - array3x1(x,y,z);

        if (numcomp > 1)
        {
            tfu.rotate(0,0, -dat2[2]);
            tfu.rotate(0, -dat2[1], 0);
            tfu.rotate(-dat2[0], 0, 0);
            
            dofu.rotate(0,0, -dat2[2]);
            dofu.rotate(0, -dat2[1], 0);
            dofu.rotate(-dat2[0], 0, 0);
        }
    }
    
    if (numcomp > 1)
    {
        tfu = tfu.resize(numcomp,1);
        dofu = dofu.resize(numcomp,1);
    }
    
    
    // Create the integration object to output:
    std::vector<integration> output(3);
    
    output[0] = integral(gamma1, dof(lambda)*tf(u));
    output[1] = integral(gamma2, -on(gamma1, invmapexpr, dof(lambda)) * factor * tfu);
    output[2] = integral(gamma1, (dof(u) - factor * on(gamma2, mapexpr, dofu)) * tf(lambda));

    return output;
}

expression sl::predefinedelasticity(expression dofu, expression tfu, expression E, expression nu, std::string myoption)
{
    // Elasticity matrix:
    expression H(6,6, {1-nu,nu,nu,0,0,0,  nu,1-nu,nu,0,0,0,  nu,nu,1-nu,0,0,0,  0,0,0,0.5*(1-2*nu),0,0,  0,0,0,0,0.5*(1-2*nu),0,  0,0,0,0,0,0.5*(1-2*nu)});
    expression coef = E/(1+nu)/(1-2*nu);
    coef.reuseit();
    H = coef * H;
    return predefinedelasticity(dofu, tfu, H, myoption);
}

expression sl::predefinedelasticity(expression dofu, expression tfu, expression H, std::string myoption)
{
    if (dofu.countrows() != tfu.countrows() || dofu.countcolumns() != 1 || tfu.countcolumns() != 1 || dofu.countrows() == 1)
    {
        std::cout << "Error in 'sl' namespace: first arguments in 'predefinedelasticity' must be either 2x1 or 3x1 vectors" << std::endl;
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
        std::cout << "Error in 'sl' namespace: invalid option or no option provided for the 2D problem in 'predefinedelasticity'" << std::endl;
        std::cout << "Available choices are: 'planestrain', 'planestress'" << std::endl;
        abort();
    }
    if (dofu.countrows() == 3)
    {
        if (myoption.length() > 0)
        {
            std::cout << "Error in 'sl' namespace: for a 3D problem the last string argument must be empty in 'predefinedelasticity'" << std::endl;
            abort();
        }
        return -( H*strain(dofu) ) * strain(tfu);
    }
    
    abort(); // fix return warning
}

expression sl::predefinedelasticity(expression dofu, expression tfu, field u, expression E, expression nu, expression prestress, std::string myoption)
{
    // Elasticity matrix:
    expression H(6,6, {1-nu,nu,nu,0,0,0,  nu,1-nu,nu,0,0,0,  nu,nu,1-nu,0,0,0,  0,0,0,0.5*(1-2*nu),0,0,  0,0,0,0,0.5*(1-2*nu),0,  0,0,0,0,0,0.5*(1-2*nu)});
    expression coef = E/(1+nu)/(1-2*nu);
    coef.reuseit();
    H = coef * H;
    return predefinedelasticity(dofu, tfu, u, H, prestress, myoption);
}

expression sl::predefinedelasticity(expression dofu, expression tfu, field u, expression H, expression prestress, std::string myoption)
{
    if (dofu.countrows() != tfu.countrows() || dofu.countcolumns() != 1 || tfu.countcolumns() != 1 || dofu.countrows() == 1)
    {
        std::cout << "Error in 'sl' namespace: first arguments in 'predefinedelasticity' must be either 2x1 or 3x1 vectors" << std::endl;
        abort();
    }
    if (dofu.countrows() == 2)
    {
        if (prestress.iszero() == false && (prestress.countcolumns() != 1 || prestress.countrows() != 3))
        {
            std::cout << "Error in 'sl' namespace: expected a 3x1 sized prestress vector (Voigt form) in 'predefinedelasticity' (set scalar 0.0 if no prestress)" << std::endl;
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
        std::cout << "Error in 'sl' namespace: invalid option or no option provided for the 2D problem in 'predefinedelasticity'" << std::endl;
        std::cout << "Available choices are: 'planestrain', 'planestress'" << std::endl;
        abort();
    }
    if (dofu.countrows() == 3)
    {
        if (prestress.iszero() == false && (prestress.countcolumns() != 1 || prestress.countrows() != 6))
        {
            std::cout << "Error in 'sl' namespace: expected a 6x1 sized prestress vector (Voigt form) in 'predefinedelasticity' (set scalar 0.0 if no prestress)" << std::endl;
            abort();
        }

        if (myoption.length() > 0)
        {
            std::cout << "Error in 'sl' namespace: for a 3D problem the last string argument must be empty in 'predefinedelasticity'" << std::endl;
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
    
    abort(); // fix return warning
}

expression sl::predefinedelectrostaticforce(expression input, expression E, expression epsilon)
{
    int md = universe::getrawmesh()->getmeshdimension();
    if (universe::isaxisymmetric)
        md++;
       
    if (md <= 1)
    {
        std::cout << "Error in 'sl' namespace: the force formula is not defined in 1D" << std::endl;
        abort();
    }
    if (input.countrows() != md)
    {
        std::cout << "Error in 'sl' namespace: the force formula expected a displacement field with " << md << " components" << std::endl;
        abort();
    }
    if (E.countrows() < md || E.countcolumns() != 1)
    {
        std::cout << "Error in 'sl' namespace: the force formula expected a " << md << "x1 E/H expression" << std::endl;
        abort();
    }
    if (epsilon.isscalar() == false)
    {
        std::cout << "Error in 'sl' namespace: the force formula is defined for a scalar epsilon/mu" << std::endl;
        abort();
    }
    E = E.resize(md,1);
    
    expression gradtfu = input;
    if (input.countcolumns() == 1)
        gradtfu = sl::grad(input);
    
    expression E2 = E*E;
        
    epsilon.reuseit();
    E.reuseit();
    E2.reuseit();
    
    expression T = epsilon * ( E*transpose(E) - 0.5*E2 * eye(md) );
    
    if (md == 2)
        T = expression(3,1,{T.at(0,0), T.at(1,1), T.at(0,1)});
    if (md == 3)
        T = expression(6,1,{T.at(0,0), T.at(1,1), T.at(2,2), T.at(1,2), T.at(0,2), T.at(0,1)});
    
    return -T*strain(gradtfu);
}

expression sl::predefinedmagnetostaticforce(expression input, expression H, expression mu)
{
    return predefinedelectrostaticforce(input, H, mu);
}

expression sl::predefinedacousticwave(expression dofp, expression tfp, expression c, expression alpha)
{
    c.reuseit(); alpha.reuseit();

    if (not(dofp.isscalar()) || not(tfp.isscalar()) || not(c.isscalar()) || not(alpha.isscalar()))
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinedacousticwave'" << std::endl;
        abort();
    }

    if (alpha.iszero())
        return ( -grad(dofp)*grad(tfp) -1.0/pow(c,2.0)*dtdt(dofp)*tfp );
    
    // Only valid for harmonic problems in case of nonzero attenuation:
    if (universe::fundamentalfrequency <= 0)
    {
        std::cout << "Error in 'sl' namespace: acoustics with nonzero attenuation is only valid for harmonic problems" << std::endl;
        abort();
    }

    return ( -grad(dofp)*grad(tfp) -1.0/pow(c,2.0)*dtdt(dofp)*tfp -2.0*alpha/c*dt(dofp)*tfp -pow(alpha,2.0)*dofp*tfp );
}

expression sl::predefinedacousticradiation(expression dofp, expression tfp, expression c, expression alpha)
{
    c.reuseit(); alpha.reuseit();

    if (not(dofp.isscalar()) || not(tfp.isscalar()) || not(c.isscalar()) || not(alpha.isscalar()))
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinedacousticradiation'" << std::endl;
        abort();
    }

    if (alpha.iszero())
        return ( -1.0/c*dt(dofp)*tfp );
    
    // Only valid for harmonic problems in case of nonzero attenuation:
    if (universe::fundamentalfrequency <= 0)
    {
        std::cout << "Error in 'sl' namespace: acoustic radiation condition with nonzero attenuation is only valid for harmonic problems" << std::endl;
        abort();
    }

    return ( -1.0/c*dt(dofp)*tfp -alpha*dofp*tfp );
}

expression sl::predefinedacousticstructureinteraction(expression dofp, expression tfp, expression dofu, expression tfu, expression c, expression rho, expression n, expression alpha, double scaling)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    
    double invscal = 1.0/scaling;

    c.reuseit(); rho.reuseit(); n.reuseit(); alpha.reuseit();

    if (not(dofp.isscalar()) || not(tfp.isscalar()) || (dofu.countcolumns() != 1 || dofu.countrows() < problemdimension) || (tfu.countcolumns() != 1 || tfu.countrows() < problemdimension) || not(c.isscalar()) || not(rho.isscalar()) || (n.countcolumns() != 1 || n.countrows() < problemdimension) || not(alpha.isscalar()))
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinedacousticstructureinteraction'" << std::endl;
        abort();
    }

    if (alpha.iszero())
        return ( -dofp*tfu*n * scaling + rho*dtdt(dofu)*n*tfp * invscal );
    
    // Only valid for harmonic problems in case of nonzero attenuation:
    if (universe::fundamentalfrequency <= 0)
    {
        std::cout << "Error in 'sl' namespace: acoustic structure interaction with nonzero attenuation is only valid for harmonic problems" << std::endl;
        abort();
    }

    return ( -dofp*tfu*n * scaling + rho*dtdt(dofu)*n*tfp * invscal +2.0*alpha*rho*c*dt(dofu)*n*tfp * invscal +rho*pow(alpha*c,2.0)*dofu*n*tfp * invscal );
}

expression sl::predefinedstokes(expression dofv, expression tfv, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant, bool isviscosityconstant)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();

    if (problemdimension < 2)
    {
        std::cout << "Error in 'sl' namespace: 'predefinedstokes' is only allowed on 2D and 3D geometries" << std::endl;
        abort();
    }
    if (dofv.countcolumns() != 1 || dofv.countrows() < problemdimension || tfv.countcolumns() != 1 || tfv.countrows() < problemdimension || not(dofp.isscalar()) || not(tfp.isscalar()) || not(mu.isscalar()) || not(rho.isscalar()))
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinedstokes'" << std::endl;
        abort();
    }

    expression output = predefinedmassconservation(dofv, tfp, rho, dtrho, gradrho, includetimederivs, isdensityconstant);

    if (includetimederivs)
        output = output - rho*dt(dofv)*tfv;

    return ( output - grad(dofp)*tfv + predefinedviscousforce(dofv, tfv, mu, isdensityconstant, isviscosityconstant) );
}

expression sl::predefinednavierstokes(expression dofv, expression tfv, expression v, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant, bool isviscosityconstant)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();

    if (problemdimension < 2)
    {
        std::cout << "Error in 'sl' namespace: 'predefinednavierstokes' is only allowed on 2D and 3D geometries" << std::endl;
        abort();
    }
    if (dofv.countcolumns() != 1 || dofv.countrows() < problemdimension || tfv.countcolumns() != 1 || tfv.countrows() < problemdimension || v.countcolumns() != 1 || v.countrows() < problemdimension || not(dofp.isscalar()) || not(tfp.isscalar()) || not(mu.isscalar()) || not(rho.isscalar()))
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinednavierstokes'" << std::endl;
        abort();
    }

    expression output = predefinedmassconservation(dofv, tfp, rho, dtrho, gradrho, includetimederivs, isdensityconstant);

    if (includetimederivs)
        output = output - rho*dt(dofv)*tfv;

    return ( output - grad(dofp)*tfv + predefinedviscousforce(dofv, tfv, mu, isdensityconstant, isviscosityconstant) +  predefinedinertialforce(dofv, tfv, v, rho) );
}


expression sl::predefinedadvectiondiffusion(expression doff, expression tff, expression v, expression alpha, expression beta, expression gamma, bool isdivvzero)
{
    int problemdimension = universe::getrawmesh()->getmeshdimension();

    bool isvsizevalid = ( v.countcolumns() == 1 && (v.iszero() || v.countrows() >= problemdimension) );

    if (not(doff.isscalar()) || not(tff.isscalar()) || not(isvsizevalid) || alpha.countrows() != alpha.countcolumns() || not(beta.isscalar()) || not(gamma.isscalar()))
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinedadvectiondiffusion'" << std::endl;
        abort();
    }

    expression output = (alpha*grad(doff)) * grad(tff);

    if (not(gamma.iszero()) && not(v.iszero()))
    {
        output = output + gamma*v*grad(doff)*tff;
        if (isdivvzero == false)
            output = output + gamma*doff*div(v)*tff;
    }

    if (beta.iszero() == false)
        output = output + beta*dt(doff)*tff;

    return output;
}

expression sl::predefineddiffusion(expression doff, expression tff, expression alpha, expression beta)
{
    return predefinedadvectiondiffusion(doff, tff, 0.0, alpha, beta, 0.0, true);
}

expression sl::predefinedstabilization(std::string stabtype, expression delta, expression f, expression v, expression diffusivity, expression residual)
{
    v.reuseit(); diffusivity.reuseit();
    
    // Avoid zero division issues for zero norm(v):
    double eps = 1e-50;
    expression invnormv = ifpositive(norm(v) - eps, 1.0/norm(v), 0.0);
    invnormv.reuseit();
    
    int problemdimension = universe::getrawmesh()->getmeshdimension();
    expression meshsize = pow(sl::meshsize(2), 1.0/problemdimension );
     
    expression doff = dof(f);
    expression tff = tf(f);
    
    if (not(residual.isscalar()) || residual.getoperationinarray(0,0)->istfincluded())
    {
        std::cout << "Error in 'sl' namespace: expected a scalar expression without test function for the residual in 'predefinedstabilization'" << std::endl;
        abort();
    }
     
    if (not(delta.isscalar()) || not(f.isscalar()) || v.countcolumns() != 1 || v.countrows() < problemdimension || diffusivity.countrows() != diffusivity.countcolumns())
    {
        std::cout << "Error in 'sl' namespace: unexpected argument dimension in 'predefinedstabilization'" << std::endl;
        abort();
    }
     
    // Isotropic diffusion term:
    if (stabtype == "iso")
    {
        expression del = delta * meshsize * norm(v);
        return ( del*grad(doff)*grad(tff) );
    }
    
    // Streamline diffusion, anisotropic:
    if (stabtype == "aniso")
    {
        expression del = delta * meshsize * invnormv;
        return ( del*(v*grad(doff))*(v*grad(tff)) );
    }
    
    // Crosswind diffusion:
    if (stabtype == "cw")
    {
        expression del = delta*pow(meshsize,1.5);
        
        expression V = eye(v.countrows()) - pow(invnormv,2.0) * v*transpose(v);
        
        return ( del * transpose(grad(doff))*V*grad(tff) );
    }
    
    // Crosswind shockwave diffusion:
    if (stabtype == "cws")
    {
        if (residual.getoperationinarray(0,0)->isdofincluded())
        {
            std::cout << "Error in 'sl' namespace: the residual cannot include a dof for cws in 'predefinedstabilization'" << std::endl;
            abort();
        }
    
        // Average diffusivity:
        expression dm = trace(diffusivity)/diffusivity.countrows();
    
        // Avoid zero division issues for zero grad(f):
        expression gradf = grad(f);
        gradf.reuseit();
        expression invnormgradf = ifpositive(norm(gradf) - eps, 1.0/norm(gradf), 0.0);
        invnormgradf.reuseit();
    
        expression vp = abs(v*gradf)*invnormgradf;
        expression gp = 0.5*meshsize*vp/dm;
        gp.reuseit();
        
        expression subtraction = ifpositive(abs(gp) - eps, delta-1.0/gp, 0.0);
        subtraction.reuseit();
        
        expression del = ifpositive(subtraction,1.0,0.0)*0.5*meshsize*subtraction*abs(residual)*invnormgradf;
        
        expression V = eye(v.countrows()) - pow(invnormv,2.0) * v*transpose(v);
        
        return ( del * transpose(grad(doff))*V*grad(tff) );
    }
    
    // Streamline diffusion, Petrov-Galerkin:
    if (stabtype == "spg")
    {
        return ( delta * meshsize * invnormv * residual*v*grad(tff) );
    }

    // Streamline diffusion, upwind Petrov-Galerkin:
    if (stabtype == "supg")
    {
        // Average diffusivity:
        expression dm = trace(diffusivity)/diffusivity.countrows();
        expression del = delta * meshsize*invnormv-dm*pow(invnormv,2.0);
        del.reuseit();
        
        expression output = del * residual*v*grad(tff);

        return ( ifpositive(del,1.0,0.0) * output );
    }

    std::cout << "Error in 'sl' namespace: unknown stabilization method '" << stabtype << "' (use 'iso', 'aniso', 'cw', 'cws', 'spg', 'supg')"  << std::endl;
    abort();
}

