#include "mathop.h"
#include "rawpoint.h"
#include "rawline.h"
#include "rawsurface.h"
#include "rawvolume.h"


double mathop::getpi(void)
{
    return 3.1415926535897932384;
}

int mathop::regionunion(const std::vector<int> physregs)
{
    universe::mymesh->getoriginalmeshpointer()->getphysicalregions()->createunion(physregs);
    return (universe::mymesh->getphysicalregions())->createunion(physregs);
}

int mathop::regionintersection(const std::vector<int> physregs)
{
    universe::mymesh->getoriginalmeshpointer()->getphysicalregions()->createintersection(physregs);
    return (universe::mymesh->getphysicalregions())->createintersection(physregs);
}

int mathop::regionall(void)
{
    universe::mymesh->getoriginalmeshpointer()->getphysicalregions()->createunionofall();
    return (universe::mymesh->getphysicalregions())->createunionofall();
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

void mathop::printvector(std::vector<bool> input)
{
    std::cout << "Vector size is " << input.size() << std::endl;
    for (int i = 0; i < input.size(); i++)
        std::cout << input[i] << " ";
    std::cout << std::endl;
}

void mathop::writevector(std::string filename, std::vector<double> towrite, char delimiter, bool writesize)
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

std::vector<double> mathop::loadvector(std::string filename, char delimiter, bool sizeincluded)
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
            myalgorithm::osclean(currentline);
            output.resize(std::stoi(currentline));
        }

        int index = 0;
        while (std::getline(name, currentline, delimiter))
        {
            myalgorithm::osclean(currentline);
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

expression mathop::norm(expression expr)
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

expression mathop::normal(int physreg)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    int elementdimension = universe::mymesh->getphysicalregions()->get(physreg)->getelementdimension();

    if (elementdimension >= problemdimension)
    {
        std::cout << "Error in 'mathop' namespace: can only compute a normal to a region of lower dimension than the geometry" << std::endl;
        abort();
    }

    expression expr;
    if (problemdimension == 1)
        return 1.0;
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

expression mathop::tangent(int physreg)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    int elementdimension = universe::mymesh->getphysicalregions()->get(physreg)->getelementdimension();

    if (elementdimension != 1 && elementdimension != 2)
    {
        std::cout << "Error in 'mathop' namespace: can only compute a tangent to a line or face region" << std::endl;
        abort();
    }

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


void mathop::setaxisymmetry(void)
{
    // Make sure the call is done before loading the mesh:
    if (universe::mymesh != NULL)
    {
        std::cout << "Error in 'mathop' namespace: 'setaxisymmetry' must be called before loading the mesh" << std::endl;
        abort();
    }
    universe::isaxisymmetric = true;
}

void mathop::setfundamentalfrequency(double f) { universe::fundamentalfrequency = f; }
void mathop::settime(double t) { universe::currenttimestep = t; }
double mathop::gettime(void) { return universe::currenttimestep; }

expression mathop::meshsize(int integrationorder)
{
    std::shared_ptr<opmeshsize> op(new opmeshsize(integrationorder));
    return expression(op);
}

expression mathop::fieldorder(field input, double alpha, double absthres)
{
    std::shared_ptr<rawfield> rf = input.getpointer();
    
    if (rf->gettypename() != "h1" && rf->gettypename() != "hcurl")
    {
        std::cout << "Error in 'mathop' namespace: field provided to 'fieldorder' must be of type 'h1' or 'hcurl' (was '" << rf->gettypename() << "')" << std::endl;
        abort();
    }
    
    if (rf->countsubfields() > 1)
    {
        std::cout << "Error in 'mathop' namespace: field provided to 'fieldorder' cannot have subfields (field of type '" << rf->gettypename(false) << "' provided has " << rf->countsubfields() << ")" << std::endl;
        std::cout << "You could instead provide the x subfield using yourfield.compx()" << std::endl;
        abort();
    }
    
    std::shared_ptr<opfieldorder> op(new opfieldorder(rf->harmonic(rf->getfirstharmonic()), alpha, absthres));

    if (alpha == -1.0)
        return expression(op);
    else
        return abs(expression(op))-1.0;
}

expression mathop::getharmonic(int harmnum, expression input, int numfftharms)
{
    if (harmnum <= 0)
    {
        std::cout << "Error in 'mathop' namespace: cannot get harmonic " << harmnum << std::endl;
        abort();
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
                std::cout << "Error in 'mathop' namespace: in 'getharmonic' expected an argument expression without dof or tf" << std::endl;
                abort();
            }

            std::shared_ptr<opharmonic> op(new opharmonic(harmnum, opin, numfftharms));
            exprs[i*n+j] = expression(op);
        }
    }

    return expression(m,n,exprs);
}

std::vector<double> mathop::gettotalforce(int physreg, expression* meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    std::vector<int> alldisjregs(universe::mymesh->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(EorH.isharmonicone(alldisjregs)) || not(epsilonormu.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'mathop' namespace: cannot have a multiharmonic argument in the total force calculation" << std::endl;
        abort();
    }
    
    int wholedomain = universe::mymesh->getphysicalregions()->createunionofall();

    int numcomps = universe::mymesh->getmeshdimension();
    if (numcomps == 1)
    {
        std::cout << "Error in 'mathop' namespace: force calculation formula is undefined in 1D" << std::endl;
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
        densematrix vecvals = fv.getvalues(intdensematrix(siz, 1, c*siz, 1));
        output[c] = vecvals.sum();
    }
    
    if (universe::isaxisymmetric)
        output = {0, 2.0*getpi()*output[1], 0};
    
    universe::mymesh->getphysicalregions()->remove({wholedomain}, false);
    
    return output;
}

std::vector<double> mathop::gettotalforce(int physreg, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    return gettotalforce(physreg, NULL, EorH, epsilonormu, extraintegrationorder);
}

std::vector<double> mathop::gettotalforce(int physreg, expression meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    return gettotalforce(physreg, &meshdeform, EorH, epsilonormu, extraintegrationorder);
}

std::vector<double> mathop::printtotalforce(int physreg, expression* meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
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

std::vector<double> mathop::printtotalforce(int physreg, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    return printtotalforce(physreg, NULL, EorH, epsilonormu, extraintegrationorder);
}

std::vector<double> mathop::printtotalforce(int physreg, expression meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder)
{
    return printtotalforce(physreg, &meshdeform, EorH, epsilonormu, extraintegrationorder);
}

void mathop::setphysicalregionshift(int shiftamount) { universe::physregshift = shiftamount; }

void mathop::writeshapefunctions(std::string filename, std::string sftypename, int elementtypenumber, int maxorder, bool allorientations)
{
    if (elementtypenumber == 7)
    {
        std::cout << "Error in 'mathop' namespace: cannot write shape functions for pyramids (non-polynomial)" << std::endl;
        abort();
    }
    if (elementtypenumber > 7)
    {
        std::cout << "Error in 'mathop' namespace: element type number must be between 0 and 7" << std::endl;
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
    std::vector<densematrix> ffvals(numcomp);

    // Prepare the element x, y and z coordinates:
    densematrix x(1,numnodes), y(1,numnodes), z(1,numnodes);
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

expression mathop::t(void) { expression exp; return exp.time(); }

void mathop::grouptimesteps(std::string filename, std::vector<std::string> filestogroup, std::vector<double> timevals)
{
    iointerface::grouptimesteps(filename, filestogroup, timevals);
}

void mathop::grouptimesteps(std::string filename, std::string fileprefix, int firstint, std::vector<double> timevals)
{
    int numsteps = timevals.size();

    std::vector<std::string> filestogroup(numsteps);
    for (int i = 0; i < numsteps; i++)
        filestogroup[i] = fileprefix+std::to_string(firstint+i)+".vtu";

    iointerface::grouptimesteps(filename, filestogroup, timevals);
}

std::vector<std::vector<shape>> mathop::loadshape(std::string meshfile)
{
    std::shared_ptr<rawmesh> loadedmesh(new rawmesh());
    
    loadedmesh->readfromfile(meshfile);
    
    nodes* loadednodes = loadedmesh->getnodes();
    std::vector<double>* nodecoords = loadednodes->getcoordinates();
    int totalnumnodes = nodecoords->size()/3;
    elements* loadedelems = loadedmesh->getelements();
    physicalregions* loadedphysregs = loadedmesh->getphysicalregions();

    int curvatureorder = loadedelems->getcurvatureorder();
    if (curvatureorder > 1)
    {
        std::cout << "Error in 'mathop' namespace: loadshape does not accept curved meshes" << std::endl;
        abort();
    }
    
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
            element myelem(j,curvatureorder);
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
            curshape = shape(std::shared_ptr<rawpoint>(new rawpoint(curphysreg, nodecoordinates, nodesinelements)));
        if (physregdim == 1)
            curshape = shape(std::shared_ptr<rawline>(new rawline(curphysreg, nodecoordinates, nodesinelements)));
        if (physregdim == 2)
            curshape = shape(std::shared_ptr<rawsurface>(new rawsurface(curphysreg, nodecoordinates, nodesinelements)));
        if (physregdim == 3)
            curshape = shape(std::shared_ptr<rawvolume>(new rawvolume(curphysreg, nodecoordinates, nodesinelements)));
    
        output[physregdim].push_back(curshape);
    }
    
    return output;
}

expression mathop::dx(expression input) { return input.spacederivative(1); }
expression mathop::dy(expression input) { return input.spacederivative(2); }
expression mathop::dz(expression input) { return input.spacederivative(3); }

expression mathop::dt(expression input) { return input.timederivative(1); }
expression mathop::dtdt(expression input) { return input.timederivative(2); }
expression mathop::dtdtdt(expression input) { return input.timederivative(3); }
expression mathop::dtdtdtdt(expression input) { return input.timederivative(4); }

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
expression mathop::exp(expression input) { return pow(2.7182818284590452353, input); }
expression mathop::mod(expression input, double modval) { return input.mod(modval); }

expression mathop::ifpositive(expression condexpr, expression trueexpr, expression falseexpr)
{
    expression output(condexpr, trueexpr, falseexpr);
    return output;
}

expression mathop::andpositive(std::vector<expression> exprs)
{
    if (exprs.size() == 0)
    {
        std::cout << "Error in 'mathop' namespace: cannot call andpositive on an empty vector of expressions" << std::endl;
        abort();
    }

    expression output(exprs[exprs.size()-1], 1, -1);

    for (int i = exprs.size()-2; i >= 0; i--)
        output = expression(exprs[i], output, -1);

    return output;
}

expression mathop::orpositive(std::vector<expression> exprs)
{
    if (exprs.size() == 0)
    {
        std::cout << "Error in 'mathop' namespace: cannot call orpositive on an empty vector of expressions" << std::endl;
        abort();
    }

    expression output(exprs[exprs.size()-1], 1, -1);

    for (int i = exprs.size()-2; i >= 0; i--)
        output = expression(exprs[i], 1, output);

    return output;
}

expression mathop::max(expression a, expression b)
{
    a.reuseit(); b.reuseit();
    
    return ifpositive(a-b, a, b);
}

expression mathop::max(field a, field b)
{
    expression expra = a;
    expression exprb = b;
    return max(expra, exprb);
}

expression mathop::max(parameter a, parameter b)
{
    expression expra = a;
    expression exprb = b;
    return max(expra, exprb);
}

expression mathop::min(expression a, expression b)
{
    a.reuseit(); b.reuseit();
    
    return ifpositive(a-b, b, a);
}

expression mathop::min(field a, field b)
{
    expression expra = a;
    expression exprb = b;
    return min(expra, exprb);
}

expression mathop::min(parameter a, parameter b)
{
    expression expra = a;
    expression exprb = b;
    return min(expra, exprb);
}


expression mathop::on(int physreg, expression expr, bool errorifnotfound) { return expr.on(physreg, NULL, errorifnotfound); }
expression mathop::on(int physreg, expression coordshift, expression expr, bool errorifnotfound) { return expr.on(physreg, &coordshift, errorifnotfound); }

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

expression mathop::eye(int size)
{
    if (size < 0)
    {
        std::cout << "Error in 'mathop' namespace: cannot create a " << size << "x" << size << " identity matrix" << std::endl;
        abort();
    }

    std::vector<expression> exprs(size);
    for (int i = 0; i < size; i++)
        exprs[i] = 1.0;
        
    expression output(size, size, exprs);
    
    return output;
}

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
        // Coordinate x on the deformed mesh:
        expression x = input.jac(2,2);
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

}

expression mathop::curl(expression input)
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
        std::cout << "Error in 'mathop' namespace: can only take the curl of an up to length 3 column vector" << std::endl;
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

expression mathop::doubledotproduct(expression a, expression b)
{
    if (a.countcolumns() != b.countcolumns() || a.countrows() != b.countrows())
    {
        std::cout << "Error in 'mathop' namespace: dimension mismatch for double dot product" << std::endl;
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

std::vector<expression> mathop::rotation(double alphax, double alphay, double alphaz, std::string type)
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
        densematrix Rx(3,3, { 1,0,0, 0,c,-s, 0,s,c });
        c = std::cos(ty); s = std::sin(ty);
        densematrix Ry(3,3, { c,0,s, 0,1,0, -s,0,c });
        c = std::cos(tz); s = std::sin(tz);
        densematrix Rz(3,3, { c,-s,0, s,c,0, 0,0,1 });
        
        densematrix R = Rz.multiply(Ry.multiply(Rx));
        
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
        densematrix Kx(6,6, { 1,0,0,0,0,0, 0,c*c,s*s,-2.0*c*s,0,0, 0,s*s,c*c,2.0*c*s,0,0, 0,c*s,-c*s,c*c-s*s,0,0, 0,0,0,0,c,s, 0,0,0,0,-s,c });
        c = std::cos(ty); s = std::sin(ty);
        densematrix Ky(6,6, { c*c,0,s*s,0,2.0*c*s,0, 0,1,0,0,0,0, s*s,0,c*c,0,-2.0*c*s,0, 0,0,0,c,0,-s, -c*s,0,c*s,0,c*c-s*s,0, 0,0,0,s,0,c });
        c = std::cos(tz); s = std::sin(tz);
        densematrix Kz(6,6, { c*c,s*s,0,0,0,-2.0*c*s, s*s,c*c,0,0,0,2.0*c*s, 0,0,1,0,0,0, 0,0,0,c,s,0, 0,0,0,-s,c,0, c*s,-c*s,0,0,0,c*c-s*s });

        densematrix K = Kz.multiply(Ky.multiply(Kx));
        
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
        densematrix invKx(6,6, { 1,0,0,0,0,0, 0,c*c,s*s,-2.0*c*s,0,0, 0,s*s,c*c,2.0*c*s,0,0, 0,c*s,-c*s,c*c-s*s,0,0, 0,0,0,0,c,s, 0,0,0,0,-s,c });
        c = std::cos(ty); s = std::sin(ty);
        densematrix invKy(6,6, { c*c,0,s*s,0,2.0*c*s,0, 0,1,0,0,0,0, s*s,0,c*c,0,-2.0*c*s,0, 0,0,0,c,0,-s, -c*s,0,c*s,0,c*c-s*s,0, 0,0,0,s,0,c });
        c = std::cos(tz); s = std::sin(tz);
        densematrix invKz(6,6, { c*c,s*s,0,0,0,-2.0*c*s, s*s,c*c,0,0,0,2.0*c*s, 0,0,1,0,0,0, 0,0,0,c,s,0, 0,0,0,-s,c,0, c*s,-c*s,0,0,0,c*c-s*s });

        densematrix invK = invKx.multiply(invKy.multiply(invKz));
        
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
    
    std::cout << "Error in 'mathop' namespace: rotation expected a type '' or 'voigt'" << std::endl;
    abort();
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

expression mathop::athp(expression expr, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt)
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

bool mathop::adapt(int verbosity)
{
    return universe::mymesh->adapthp(verbosity);
}

expression mathop::zienkiewiczzhu(expression input)
{
    std::vector<int> alldisjregs(universe::mymesh->getdisjointregions()->count());
    std::iota(alldisjregs.begin(), alldisjregs.end(), 0);
    if (not(input.isharmonicone(alldisjregs)))
    {
        std::cout << "Error in 'mathop' namespace: in 'zienkiewiczzhu' cannot have a multiharmonic expression as argument" << std::endl;
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

    vec sol(std::shared_ptr<rawvec>(new rawvec(b.getpointer()->getdofmanager())));
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

std::vector<vec> mathop::solve(mat A, std::vector<vec> b, std::string soltype)
{
    if (soltype != "lu")
    {
        std::cout << "Error in 'mathop' namespace: unknown direct solver type '" << soltype << "' (use 'lu')" << std::endl;
        abort();
    }
    for (int i = 0; i < b.size(); i++)
    {
        if (A.countrows() != b[i].size())
        {
            std::cout << "Error in 'mathop' namespace: multi-rhs direct solve of Ax = b failed (size of A and at least one rhs do not match)" << std::endl;
            abort();
        }
        if (A.getpointer() == NULL || b[i].getpointer() == NULL)
        {
            std::cout << "Error in 'mathop' namespace: multi-rhs direct solve of Ax = b failed (A or at least one rhs is undefined)" << std::endl;
            abort();
        }
    }
    
    if (b.size() == 0)
        return {};
        
    int numrhs = b.size();
    int len = b[0].size();
    
    // The copy of the rhs is returned in case there is no nonzero entry in A:
    if (A.countnnz() == 0)
    {
        std::vector<vec> bcopy(numrhs);
        for (int i = 0; i < numrhs; i++)
            bcopy[i] = b[i].copy();
        return bcopy;
    }

    // Concatenate rhs vecs to densematrix:
    intdensematrix ads(len, 1, 0, 1);
    densematrix rhs(numrhs, len);
    double* rhsptr = rhs.getvalues();
    for (int i = 0; i < numrhs; i++)
    {
        densematrix vecvals = b[i].getvalues(ads);
        double* vecvalsptr = vecvals.getvalues();
        for (int j = 0; j < len; j++)
            rhsptr[i*len+j] = vecvalsptr[j];
    }
    
    // Solve multi-rhs:
    densematrix sols = solve(A, rhs, soltype);
    double* solsptr = sols.getvalues();

    // Extract 'sols' rows to sol vecs:
    densematrix vals(len,1);
    double* valsptr = vals.getvalues();
    std::vector<vec> outvecs(numrhs);
    for (int i = 0; i < numrhs; i++)
    {
        for (int j = 0; j < len; j++)
            valsptr[j] = solsptr[i*len+j];
    
        outvecs[i] = vec(std::shared_ptr<rawvec>(new rawvec(b[i].getpointer()->getdofmanager())));
        outvecs[i].setvalues(ads, vals);
    }
    
    return outvecs;
}

densematrix mathop::solve(mat A, densematrix b, std::string soltype)
{
    int numrhs = b.countrows();
    int len = b.countcolumns();
 
    Mat Apetsc = A.getpetsc();
    
    KSP* ksp = A.getpointer()->getksp();
    PC pc;
    if (A.getpointer()->isludefined() == false)
    {
        KSPCreate(PETSC_COMM_WORLD, ksp);
        KSPSetOperators(*ksp, Apetsc, Apetsc);
        KSPSetFromOptions(*ksp);

        KSPGetPC(*ksp,&pc);
        PCSetType(pc,PCLU);
        PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
        PCSetUp(pc);
    }
    else
        KSPGetPC(*ksp,&pc);
        
    PCFactorGetMatrix(pc, &Apetsc);

    densematrix densesols(numrhs, len);

    Mat sols, rhses;
    MatCreateSeqDense(PETSC_COMM_SELF, len, numrhs, densesols.getvalues(), &sols);
    MatCreateSeqDense(PETSC_COMM_SELF, len, numrhs, b.getvalues(), &rhses);
    
    // 'rhs' and 'sols' are considered column major in petsc.
    MatMatSolve(Apetsc, rhses, sols);
    
    MatDestroy(&sols);
    MatDestroy(&rhses);

    A.getpointer()->isludefined(true);

    if (A.getpointer()->islutobereused() == false)
    {
        KSPDestroy(ksp);
        A.getpointer()->isludefined(false);
    }
    
    return densesols;
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
    if (precondtype != "ilu" && precondtype != "sor" && precondtype != "gamg")
    {
        std::cout << "Error in 'mathop' namespace: unknown preconditioner type '" << precondtype << "' (use 'ilu', 'sor' or 'gamg')" << std::endl;
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
    if (precondtype == "gamg")
        PCSetType(pc,PCGAMG);

    KSPSolve(*ksp, bpetsc, solpetsc);

    // Get the number of required iterations and the residual norm:
    KSPGetIterationNumber(*ksp, &maxnumit);
    KSPGetResidualNorm(*ksp, &relrestol);

    KSPDestroy(ksp);
}

void mathop::solve(formulation formul)
{
    // Remove leftovers (if any):
    mat A = formul.A(); vec b = formul.b();
    // Generate:
    formul.generate();
    // Solve:
    vec sol = mathop::solve(formul.A(), formul.b());

    // Save to fields:
    setdata(sol);
}

void mathop::solve(std::vector<formulation> formuls)
{
    for (int i = 0; i < formuls.size(); i++)
        solve(formuls[i]);
}


std::vector<double> mathop::linspace(double a, double b, int num)
{
    if (num < 0)
    {
        std::cout << "Error in 'mathop' namespace: cannot call 'linspace' for " << num << " points" << std::endl;
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

std::vector<double> mathop::logspace(double a, double b, int num, double basis)
{
    if (num < 0)
    {
        std::cout << "Error in 'mathop' namespace: cannot call 'logspace' for " << num << " points" << std::endl;
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

expression mathop::dbtoneper(expression toconvert)
{
    if (toconvert.iszero())
        return toconvert;

    return ( 0.1151292546497023 * toconvert );
}

void mathop::setdata(vec invec)
{
    if (invec.getpointer() == NULL)
        return;
    
    // Get all fields in the vec structure:
    std::vector<std::shared_ptr<rawfield>> allfields = invec.getpointer()->getdofmanager()->getfields();
    
    for (int i = 0; i < allfields.size(); i++)
        allfields[i]->setdata(-1, invec|field(allfields[i]));
}


////////// PREDEFINED OPERATORS

expression mathop::strain(expression input)
{
    if ((input.countrows() != 2 && input.countrows() != 3) || (input.countcolumns() != 1 && input.countrows() != input.countcolumns()))
    {
        std::cout << "Error in 'mathop' namespace: can only compute the strains of a 2x1 or 3x1 column vector or its gradient" << std::endl;
        abort();
    }

    expression gradu = input;
    if (input.countcolumns() == 1)
        gradu = mathop::grad(input);

    if (input.countrows() == 2)
        return expression(3,1,{gradu.at(0,0), gradu.at(1,1), gradu.at(1,0) + gradu.at(0,1)});
    if (input.countrows() == 3)
        return expression(6,1,{gradu.at(0,0), gradu.at(1,1), gradu.at(2,2), gradu.at(2,1) + gradu.at(1,2), gradu.at(0,2) + gradu.at(2,0), gradu.at(0,1) + gradu.at(1,0)});
}

expression mathop::greenlagrangestrain(expression input)
{
    if ((input.countrows() != 2 && input.countrows() != 3) || (input.countcolumns() != 1 && input.countrows() != input.countcolumns()))
    {
        std::cout << "Error in 'mathop' namespace: can only compute the green-lagrange strains of a 2x1 or 3x1 column vector or its gradient" << std::endl;
        abort();
    }

    expression gradu = input;
    if (input.countcolumns() == 1)
        gradu = mathop::grad(input);

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
        return ( - mu*doubledotproduct(grad(dofv), grad(tfv)) );

    if (isdensityconstant)
        return ( - mu*doubledotproduct(grad(dofv), grad(tfv)) - mu*doubledotproduct(transpose(grad(dofv)), grad(tfv)) );
    else
        return ( - mu*doubledotproduct(grad(dofv), grad(tfv)) - mu*doubledotproduct(transpose(grad(dofv)), grad(tfv)) + (2.0/3.0)*mu*div(dofv)*trace(grad(tfv)) );
}


////////// PREDEFINED FORMULATIONS

std::vector<integration> mathop::continuitycondition(int gamma1, int gamma2, field u1, field u2, int lagmultorder, bool errorifnotfound)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    int gamma1dim = universe::mymesh->getphysicalregions()->get(gamma1)->getelementdimension();
    int gamma2dim = universe::mymesh->getphysicalregions()->get(gamma2)->getelementdimension();
    if (gamma1dim != gamma2dim || gamma1dim >= problemdimension)
    {
        std::cout << "Error in 'mathop' namespace: expected boundary regions for gamma1 and gamma2 in 'continuitycondition'" << std::endl;
        abort();
    }
    
    std::shared_ptr<rawfield> ptr1 = u1.getpointer();
    std::shared_ptr<rawfield> ptr2 = u2.getpointer();
    
    // Make sure the fields are similar:
    if (ptr1->gettypename(false) != ptr2->gettypename(false) || ptr1->getharmonics() != ptr2->getharmonics())
    {
        std::cout << "Error in 'mathop' namespace: in 'continuitycondition' expected two fields of same type and harmonic content" << std::endl;
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

std::vector<integration> mathop::continuitycondition(int gamma1, int gamma2, field u1, field u2, std::vector<double> rotcent, double rotangz, double angzmod, double factor, int lagmultorder)
{       
    int problemdimension = universe::mymesh->getmeshdimension();
    int gamma1dim = universe::mymesh->getphysicalregions()->get(gamma1)->getelementdimension();
    int gamma2dim = universe::mymesh->getphysicalregions()->get(gamma2)->getelementdimension();
    if (gamma1dim != gamma2dim || gamma1dim >= problemdimension)
    {
        std::cout << "Error in 'mathop' namespace: expected boundary regions for gamma1 and gamma2 in 'continuitycondition'" << std::endl;
        abort();
    }

    std::shared_ptr<rawfield> ptr1 = u1.getpointer();
    std::shared_ptr<rawfield> ptr2 = u2.getpointer();
    
    // Make sure the fields are similar:
    if (ptr1->gettypename(false) != ptr2->gettypename(false) || ptr1->getharmonics() != ptr2->getharmonics())
    {
        std::cout << "Error in 'mathop' namespace: in 'continuitycondition' expected two fields of same type and harmonic content" << std::endl;
        abort();
    }
    
    if (rotcent.size() != 3)
    {
        std::cout << "Error in 'mathop' namespace: in 'continuitycondition' expected a vector of length 3 as fifth argument" << std::endl;
        abort();
    }
    
    if (factor != -1 && factor != 1)
    {
        std::cout << "Error in 'mathop' namespace: in 'continuitycondition' the factor must be -1 or 1" << std::endl;
        abort();
    }

    if (angzmod < 0.0 || angzmod > 180.0)
    {
        std::cout << "Error in 'mathop' namespace: in 'continuitycondition' the angular modulo should be in range [0,180]" << std::endl;
        abort();
    }

    if (rotangz < 0.0 || rotangz > angzmod)
    {
        std::cout << "Error in 'mathop' namespace: in 'continuitycondition' the rotation angle should be in range [0,angzmod]" << std::endl;
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

std::vector<integration> mathop::periodicitycondition(int gamma1, int gamma2, field u, std::vector<double> dat1, std::vector<double> dat2, double factor, int lagmultorder)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    int gamma1dim = universe::mymesh->getphysicalregions()->get(gamma1)->getelementdimension();
    int gamma2dim = universe::mymesh->getphysicalregions()->get(gamma2)->getelementdimension();
    if (gamma1dim != gamma2dim || gamma1dim >= problemdimension)
    {
        std::cout << "Error in 'mathop' namespace: expected boundary regions for gamma1 and gamma2 in 'periodicitycondition'" << std::endl;
        abort();
    }
    
    if (dat1.size() != 3)
    {
        std::cout << "Error in 'mathop' namespace: in 'periodicitycondition' expected a vector of length 3 as fourth argument" << std::endl;
        abort();
    }
    if (dat2.size() != 1 && dat2.size() != 3)
    {
        std::cout << "Error in 'mathop' namespace: in 'periodicitycondition' expected a vector of length 1 or 3 as fifth argument" << std::endl;
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

expression mathop::predefinedelasticity(expression dofu, expression tfu, expression E, expression nu, std::string myoption)
{
    // Elasticity matrix:
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
    // Elasticity matrix:
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

expression mathop::predefinedacoustics(expression dofp, expression tfp, expression c, expression alpha)
{
    c.reuseit(); alpha.reuseit();

    if (not(dofp.isscalar()) || not(tfp.isscalar()) || not(c.isscalar()) || not(alpha.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinedacoustics'" << std::endl;
        abort();
    }

    if (alpha.iszero())
        return ( -grad(dofp)*grad(tfp) -1.0/pow(c,2.0)*dtdt(dofp)*tfp );
    
    // Only valid for harmonic problems in case of nonzero attenuation:
    if (universe::fundamentalfrequency == -1)
    {
        std::cout << "Error in 'mathop' namespace: acoustics with nonzero attenuation is only valid for harmonic problems" << std::endl;
        abort();
    }

    return ( -grad(dofp)*grad(tfp) -1.0/pow(c,2.0)*dtdt(dofp)*tfp -2.0*alpha/c*dt(dofp)*tfp -pow(alpha,2.0)*dofp*tfp );
}

expression mathop::predefinedacousticradiation(expression dofp, expression tfp, expression c, expression alpha)
{
    c.reuseit(); alpha.reuseit();

    if (not(dofp.isscalar()) || not(tfp.isscalar()) || not(c.isscalar()) || not(alpha.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinedacousticradiation'" << std::endl;
        abort();
    }

    if (alpha.iszero())
        return ( -1.0/c*dt(dofp)*tfp );
    
    // Only valid for harmonic problems in case of nonzero attenuation:
    if (universe::fundamentalfrequency == -1)
    {
        std::cout << "Error in 'mathop' namespace: acoustic radiation condition with nonzero attenuation is only valid for harmonic problems" << std::endl;
        abort();
    }

    return ( -1.0/c*dt(dofp)*tfp -alpha*dofp*tfp );
}

expression mathop::predefinedacousticstructureinteraction(expression dofp, expression tfp, expression dofu, expression tfu, expression c, expression rho, expression n, expression alpha, double scaling)
{
    int problemdimension = universe::mymesh->getmeshdimension();
    
    double invscal = 1.0/scaling;

    c.reuseit(); rho.reuseit(); n.reuseit(); alpha.reuseit();

    if (not(dofp.isscalar()) || not(tfp.isscalar()) || (dofu.countcolumns() != 1 || dofu.countrows() < problemdimension) || (tfu.countcolumns() != 1 || tfu.countrows() < problemdimension) || not(c.isscalar()) || not(rho.isscalar()) || (n.countcolumns() != 1 || n.countrows() < problemdimension) || not(alpha.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinedacousticstructureinteraction'" << std::endl;
        abort();
    }

    if (alpha.iszero())
        return ( -dofp*tfu*n * scaling + rho*dtdt(dofu)*n*tfp * invscal );
    
    // Only valid for harmonic problems in case of nonzero attenuation:
    if (universe::fundamentalfrequency == -1)
    {
        std::cout << "Error in 'mathop' namespace: acoustic structure interaction with nonzero attenuation is only valid for harmonic problems" << std::endl;
        abort();
    }

    return ( -dofp*tfu*n * scaling + rho*dtdt(dofu)*n*tfp * invscal +2.0*alpha*rho*c*dt(dofu)*n*tfp * invscal +rho*pow(alpha*c,2.0)*dofu*n*tfp * invscal );
}

expression mathop::predefinedstokes(expression dofv, expression tfv, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant, bool isviscosityconstant)
{
    int problemdimension = universe::mymesh->getmeshdimension();

    if (problemdimension < 2)
    {
        std::cout << "Error in 'mathop' namespace: 'predefinedstokes' is only allowed on 2D and 3D geometries" << std::endl;
        abort();
    }
    if (dofv.countcolumns() != 1 || dofv.countrows() < problemdimension || tfv.countcolumns() != 1 || tfv.countrows() < problemdimension || not(dofp.isscalar()) || not(tfp.isscalar()) || not(mu.isscalar()) || not(rho.isscalar()))
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
    if (dofv.countcolumns() != 1 || dofv.countrows() < problemdimension || tfv.countcolumns() != 1 || tfv.countrows() < problemdimension || v.countcolumns() != 1 || v.countrows() < problemdimension || not(dofp.isscalar()) || not(tfp.isscalar()) || not(mu.isscalar()) || not(rho.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinednavierstokes'" << std::endl;
        abort();
    }

    expression output = predefinedmassconservation(dofv, tfp, rho, dtrho, gradrho, includetimederivs, isdensityconstant);

    if (includetimederivs)
        output = output - rho*dt(dofv)*tfv;

    return ( output - grad(dofp)*tfv + predefinedviscousforce(dofv, tfv, mu, isdensityconstant, isviscosityconstant) +  predefinedinertialforce(dofv, tfv, v, rho) );
}


expression mathop::predefinedadvectiondiffusion(expression doff, expression tff, expression v, expression alpha, expression beta, expression gamma, bool isdivvzero)
{
    int problemdimension = universe::mymesh->getmeshdimension();

    bool isvsizevalid = ( v.countcolumns() == 1 && (v.iszero() || v.countrows() >= problemdimension) );

    if (not(doff.isscalar()) || not(tff.isscalar()) || not(isvsizevalid) || alpha.countrows() != alpha.countcolumns() || not(beta.isscalar()) || not(gamma.isscalar()))
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinedadvectiondiffusion'" << std::endl;
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

expression mathop::predefineddiffusion(expression doff, expression tff, expression alpha, expression beta)
{
    return predefinedadvectiondiffusion(doff, tff, 0.0, alpha, beta, 0.0, true);
}

expression mathop::predefinedstabilization(std::string stabtype, expression delta, expression f, expression v, expression diffusivity, expression residual)
{
    v.reuseit(); diffusivity.reuseit();
    
    // Avoid zero division issues for zero norm(v):
    double eps = 1e-50;
    expression invnormv = ifpositive(norm(v) - eps, 1.0/norm(v), 0.0);
    invnormv.reuseit();
    
    int problemdimension = universe::mymesh->getmeshdimension();
    expression meshsize = pow(mathop::meshsize(2), 1.0/problemdimension );
     
    expression doff = dof(f);
    expression tff = tf(f);
    
    if (not(residual.isscalar()) || residual.getoperationinarray(0,0)->istfincluded())
    {
        std::cout << "Error in 'mathop' namespace: expected a scalar expression without test function for the residual in 'predefinedstabilization'" << std::endl;
        abort();
    }
     
    if (not(delta.isscalar()) || not(f.isscalar()) || v.countcolumns() != 1 || v.countrows() < problemdimension || diffusivity.countrows() != diffusivity.countcolumns())
    {
        std::cout << "Error in 'mathop' namespace: unexpected argument dimension in 'predefinedstabilization'" << std::endl;
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
            std::cout << "Error in 'mathop' namespace: the residual cannot include a dof for cws in 'predefinedstabilization'" << std::endl;
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

    std::cout << "Error in 'mathop' namespace: unknown stabilization method '" << stabtype << "' (use 'iso', 'aniso', 'cw', 'cws', 'spg', 'supg')"  << std::endl;
    abort();
}

