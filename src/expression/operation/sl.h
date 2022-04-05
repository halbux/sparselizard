// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.
//
// Thanks to R. Haouari for the stabilization terms.

#ifndef SL_H
#define SL_H

#include <iostream>
#include <string>
#include "expression.h"
#include "integration.h"
#include "universe.h"
#include "iointerface.h"
#include "vec.h"
#include "mat.h"
#include "formulation.h"
#include "rawmesh.h"
#include "dofmanager.h"

class rawmesh;
class expression;
class integration;
class mat;
class vec;
class field;
class parameter;
class formulation;
class shape;
class dofmanager;

namespace sl
{
    // Get the code version number and name:
    void printversion(void);
    int getversion(void);
    int getsubversion(void);
    std::string getversionname(void);
    
    void setmaxnumthreads(int mnt);
    int getmaxnumthreads(void);
    
    double getpi(void);

    // Get a random value uniformly distributed between 0.0 and 1.0:
    double getrandom(void);
    
    // Perform operations (union, intersection) on physical regions:
    int selectunion(std::vector<int> physregs);
    int selectintersection(std::vector<int> physregs);
    int selectall(void);
    
    // Check if a region is defined/empty/fully included in another region/touches another region:
    bool isdefined(int physreg);
    bool isempty(int physreg);
    bool isinside(int physregtocheck, int physreg);
    bool istouching(int physregtocheck, int physreg);

    void printvector(std::vector<double> input);
    void printvector(std::vector<int> input);
    void printvector(std::vector<bool> input);

    void writevector(std::string filename, std::vector<double> towrite, char delimiter = ',', bool writesize = false);
    // Load a vector of doubles separated by a character:
    std::vector<double> loadvector(std::string filename, char delimiter = ',', bool sizeincluded = false);

    // Partition the mesh into numranks parts and return the mesh file name for each rank:
    std::string allpartition(std::string meshfile);

    // Compute the L2 norm of an expression:
    expression norm(expression expr);

    // Normal vector with unit norm and pointing outward of a physical region:
    expression normal(void);
    expression normal(int physreg);
    expression getnormal(int physreg);
    // Tangent vector with unit norm:
    expression tangent(void);

    // Write scalar or vector values at given coordinates to file:
    void scatterwrite(std::string filename, std::vector<double> xcoords, std::vector<double> ycoords, std::vector<double> zcoords, std::vector<double> compxevals, std::vector<double> compyevals = {}, std::vector<double> compzevals = {});

    void setaxisymmetry(void);

    void setfundamentalfrequency(double f);
    void settime(double t);
    double gettime(void);
    
    expression meshsize(int integrationorder);
    // Return the field order to hold alpha % of the total coefficient weight. Return the actual field order with alpha set to -1.0.
    expression fieldorder(field input, double alpha = -1.0, double absthres = 0.0);
    
    // Get a single harmonic:
    expression getharmonic(int harmnum, expression input, int numfftharms = -1);
    // Make a harmonic expression:
    expression makeharmonic(std::vector<int> harms, std::vector<expression> exprs);
    // Return an expression containing the origin harmonics moved to new harmonic positions:
    expression moveharmonic(std::vector<int> origharms, std::vector<int> destharms, expression input, int numfftharms = -1);
    
    std::vector<double> gettotalforce(int physreg, expression* meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder);
    std::vector<double> gettotalforce(int physreg, expression EorH, expression epsilonormu, int extraintegrationorder = 0);
    std::vector<double> gettotalforce(int physreg, expression meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder = 0);
    
    std::vector<double> printtotalforce(int physreg, expression* meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder);
    std::vector<double> printtotalforce(int physreg, expression EorH, expression epsilonormu, int extraintegrationorder = 0);
    std::vector<double> printtotalforce(int physreg, expression meshdeform, expression EorH, expression epsilonormu, int extraintegrationorder = 0);
    
    void setphysicalregionshift(int shiftamount);

    // Write all shape functions for an element type up to a given order:
    void writeshapefunctions(std::string filename, std::string sftypename, int elementtypenumber, int maxorder, bool allorientations = false);
    
    // The time variable:
    expression t(void);

    // Group .vtu timestep files in a .pvd file:
    void grouptimesteps(std::string filename, std::vector<std::string> filestogroup, std::vector<double> timevals);
    void grouptimesteps(std::string filename, std::string fileprefix, int firstint, std::vector<double> timevals);
    
    // Load a mesh file into a shape. Return the list of shapes of all dimensions {{0D},{1D},{2D},{3D}} (if any):
    std::vector<std::vector<shape>> loadshape(std::string meshfile);
    
    // Make the time derivative vectors available in the universe:
    void settimederivative(vec dtx);
    void settimederivative(vec dtx, vec dtdtx);

    expression dx(expression input);
    expression dy(expression input);
    expression dz(expression input);

    expression dt(expression input);
    expression dtdt(expression input);
    expression dtdtdt(expression input);
    expression dtdtdtdt(expression input);
    
    // Transient approximation of dt and dtdt:
    expression dt(expression input, double initdt, double initdtdt);
    expression dtdt(expression input, double initdt, double initdtdt);
    
    expression transientdtapprox(int dtorder, expression input, double initdt, double initdtdt);

    expression sin(expression input);
    expression cos(expression input);
    expression tan(expression input);
    expression asin(expression input);
    expression acos(expression input);
    expression atan(expression input);
    expression abs(expression input);
    expression sqrt(expression input);
    expression log10(expression input);
    expression pow(expression base, expression exponent);
    expression exp(expression input);
    expression mod(expression input, double modval);

    // Easy conditional functions for expressions (true if the expression value is >= 0):
    expression ifpositive(expression condexpr, expression trueexpr, expression falseexpr);
    expression andpositive(std::vector<expression> exprs);
    expression orpositive(std::vector<expression> exprs);
    
    expression max(expression a, expression b);
    expression max(field a, field b);
    expression max(parameter a, parameter b);
    expression min(expression a, expression b);
    expression min(field a, field b);
    expression min(parameter a, parameter b);
    

    // Evaluate an expression on physical region 'physreg' using interpolation:
    expression on(int physreg, expression expr, bool errorifnotfound = true);
    // Interpolate at coordinates shifted by 'coordshift':
    expression on(int physreg, expression coordshift, expression expr, bool errorifnotfound = true);

    expression comp(int selectedcomp, expression input);
    expression compx(expression input);
    expression compy(expression input);
    expression compz(expression input);

    expression entry(int row, int col, expression input);
    
    expression eye(int size);

    expression transpose(expression input);
    expression inverse(expression input);
    expression determinant(expression input);

    expression grad(expression input);
    expression div(expression input);
    expression curl(expression input);

    // Cross product between vector a and vector b.
    // Any argument that does not have 3 components will be filled with zeros.
    expression crossproduct(expression a, expression b);

    // Double dot a:b product:
    expression doubledotproduct(expression a, expression b);

    // Get the trace of a square matrix:
    expression trace(expression a);
    
    // Get the rotation matrix and its inverse for a 3x3 tensor with (x,y,z) component 
    // ordering or a 6x6 tensor in Voigt form. Input angles are in degrees.
    std::vector<expression> rotation(double alphax, double alphay, double alphaz, std::string type = "");

    integration integral(int physreg, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    integration integral(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    // For the multiharmonic resolution an extra integer is required
    // to know with how many harmonics the coef multiplying the tf
    // and/or dof should be approximated. Set 'numcoefharms' negative
    // and it will be as if you were calling the above non
    // multiharmonic functions.
    integration integral(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    integration integral(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);

    expression dof(expression input);
    expression dof(expression input, int physreg);
    expression tf(expression input);
    expression tf(expression input, int physreg);
    
    // Return an expression whose value will be calculated on the argument mesh state (the expression is hp-synchronized to that mesh state after the call):
    expression athp(expression expr, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt);

    // hp-adaptation:
    bool adapt(int verbosity = 0);
    bool alladapt(int verbosity = 0);
    
    // Define a Zienkiewicz-Zhu type error indicator:
    expression zienkiewiczzhu(expression input);

    // Define typically used arrays for convenience:
    expression array1x1(expression term11);
    expression array1x2(expression term11, expression term12);
    expression array1x3(expression term11, expression term12, expression term13);
    expression array2x1(expression term11, expression term21);
    expression array2x2(expression term11, expression term12, expression term21, expression term22);
    expression array2x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23);
    expression array3x1(expression term11, expression term21, expression term31);
    expression array3x2(expression term11, expression term12, expression term21, expression term22, expression term31, expression term32);
    expression array3x3(expression term11, expression term12, expression term13, expression term21, expression term22, expression term23, expression term31, expression term32, expression term33);


    // Direct resolution (with or without diagonal scaling):
    vec solve(mat A, vec b, std::string soltype = "lu", bool diagscaling = false);
    // Multi-rhs direct resolution:
    std::vector<vec> solve(mat A, std::vector<vec> b, std::string soltype = "lu");
    // Densematrix 'b' has size #rhs x #dofs:
    densemat solve(mat A, densemat b, std::string soltype);
    
    // Iterative resolution (with or without diagonal scaling):
    void solve(mat A, vec b, vec sol, double& relrestol, int& maxnumit, std::string soltype = "bicgstab", std::string precondtype = "sor", int verbosity = 1, bool diagscaling = false);
    
    
    // Exchange densemat data with MPI:
    void exchange(std::vector<int> targetranks, std::vector<densemat> sends, std::vector<densemat> receives);
    
    // MPI based gmres with custom matrix free product (no restart). Initial guess and solution are in x.
    // Relative residual at each iteration is returned. Length is number of iterations + 1 (first is initial residual).
    std::vector<double> gmres(densemat (*mymatmult)(densemat), densemat b, densemat x, double relrestol, int maxnumit, int verbosity = 1);
    
    // Know which dofs to send and at which dofs to receive for the DDM. Choose the rawfields and the domain interface dimensions (length 3) to consider.
    void mapdofs(std::shared_ptr<dofmanager> dm, std::vector<std::shared_ptr<rawfield>> rfs, std::vector<bool> isdimactive, std::vector<indexmat>& sendinds, std::vector<indexmat>& recvinds);
    
    std::vector<double> linspace(double a, double b, int num);
    std::vector<double> logspace(double a, double b, int num, double basis = 10.0);
    
    // Convert dB to Nepers:
    expression dbtoneper(expression toconvert);
    
    
    // Set all fields and ports to the values available in the vec object:
    void setdata(vec invec);


    ////////// PREDEFINED OPERATORS

    // Gives the engineering strains of a 2D or 3D mechanical displacement vector (Voigt form):
    expression strain(expression input);
    // Gives the Green-Lagrange strains of a 2D or 3D mechanical displacement vector (Voigt form):
    expression greenlagrangestrain(expression input);
    // Gives the von Mises stress (Voigt form expected for the argument):
    expression vonmises(expression stress);

    // Weak form of the mass conservation for Navier-Stokes:
    expression predefinedmassconservation(expression dofv, expression tfp, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant);
    // Weak form of the inertial forces for Navier-Stokes:
    expression predefinedinertialforce(expression dofv, expression tfv, expression v, expression rho);
    // Weak form of the viscous forces for Navier-Stokes:
    expression predefinedviscousforce(expression dofv, expression tfv, expression mu, bool isdensityconstant, bool isviscosityconstant);


    ////////// PREDEFINED FORMULATIONS
    
    std::vector<integration> continuitycondition(int gamma1, int gamma2, field u1, field u2, int lagmultorder, bool errorifnotfound = true);
    std::vector<integration> continuitycondition(int gamma1, int gamma2, field u1, field u2, std::vector<double> rotcent, double rotangz, double angzmod, double factor, int lagmultorder);
    std::vector<integration> periodicitycondition(int gamma1, int gamma2, field u, std::vector<double> dat1, std::vector<double> dat2, double factor, int lagmultorder);

    // Isotropic linear elasticity:
    expression predefinedelasticity(expression dofu, expression tfu, expression Eyoung, expression nupoisson, std::string myoption = "");
    // General anisotropic linear elasticity:
    expression predefinedelasticity(expression dofu, expression tfu, expression elasticitymatrix, std::string myoption = "");

    // Isotropic elasticity with geometrical nonlinearity and prestress (ignored if zero):
    expression predefinedelasticity(expression dofu, expression tfu, field u, expression Eyoung, expression nupoisson, expression prestress, std::string myoption = "");
    // General anisotropic elasticity with geometrical nonlinearity and prestress (ignored if zero):
    expression predefinedelasticity(expression dofu, expression tfu, field u, expression elasticitymatrix, expression prestress, std::string myoption = "");

    expression predefinedelectrostaticforce(expression input, expression E, expression epsilon);
    expression predefinedmagnetostaticforce(expression input, expression H, expression mu);

    expression predefinedacousticwave(expression dofp, expression tfp, expression soundspeed, expression neperattenuation);
    expression predefinedacousticradiation(expression dofp, expression tfp, expression soundspeed, expression neperattenuation);
    expression predefinedacousticstructureinteraction(expression dofp, expression tfp, expression dofu, expression tfu, expression soundspeed, expression fluiddensity, expression normal, expression neperattenuation, double scaling = 1.0);

    // Stokes flow for Newtonian fluids:
    expression predefinedstokes(expression dofv, expression tfv, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs = false, bool isdensityconstant = true, bool isviscosityconstant = true);
    // Navier-Stokes flow for Newtonian fluids:
    expression predefinednavierstokes(expression dofv, expression tfv, expression v, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs = false, bool isdensityconstant = true, bool isviscosityconstant = true);
    
    // Advection-diffusion equations:
    expression predefinedadvectiondiffusion(expression doff, expression tff, expression v, expression alpha, expression beta, expression gamma, bool isdivvzero = true);
    expression predefineddiffusion(expression doff, expression tff, expression alpha, expression beta);
    
    // Stabilization for advection-diffusion problems:
    expression predefinedstabilization(std::string stabtype, expression delta, expression f, expression v, expression diffusivity, expression residual);
};

#endif
