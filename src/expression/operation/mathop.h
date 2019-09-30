// sparselizard - Copyright (C) 2017- A. Halbach
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef MATHOP_H
#define MATHOP_H

#include <iostream>
#include <string>
#include "expression.h"
#include "integration.h"
#include "universe.h"
#include "iointerface.h"
#include "vec.h"
#include "mat.h"
#include "formulation.h"

class expression;
class integration;
class mat;
class vec;
class field;
class formulation;
class shape;

namespace mathop
{
    double getpi(void);

    // Perform operations (union, intersection, exclusion) on physical regions:
    int regionunion(const std::vector<int> physregs);
    int regionintersection(const std::vector<int> physregs);
    int regionexclusion(int physreg, int toexclude);

    void printvector(std::vector<double> input);
    void printvector(std::vector<int> input);

    void writevector(std::string filename, std::vector<double> towrite, char delimiter = ',', bool writesize = false);
    // Load a vector of doubles separated by a character:
    std::vector<double> loadvector(std::string filename, char delimiter = ',', bool sizeincluded = false);

    // Compute the L2 norm of a vector expression:
    expression norm(expression expr);

    // Define the (unit norm) vector normal to a surface in 3D (or to a line in 2D).
    expression normal(int surfphysreg);

    // Write scalar or vector values at given coordinates to file:
    void scatterwrite(std::string filename, std::vector<double> xcoords, std::vector<double> ycoords, std::vector<double> zcoords, std::vector<double> compxevals, std::vector<double> compyevals = {}, std::vector<double> compzevals = {});

    void setaxisymmetry(void);

    void setfundamentalfrequency(double f);
    void settime(double t);
    double gettime(void);

    // The time variable:
    expression t(void);

    // Group .vtu timestep files in a .pvd file:
    void grouptimesteps(std::string filename, std::vector<std::string> filestogroup, std::vector<double> timevals);
    void grouptimesteps(std::string filename, std::string fileprefix, int firstint, std::vector<double> timevals);
    
    // Load a mesh file into a shape. Return the list of shapes of all dimensions {{0D},{1D},{2D},{3D}} (if any):
    std::vector<std::vector<shape>> loadshape(std::string meshfile);

    expression dx(expression input);
    expression dy(expression input);
    expression dz(expression input);

    expression dt(expression input);
    expression dtdt(expression input);
    expression dtdtdt(expression input);
    expression dtdtdtdt(expression input);

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

    // Easy conditional functions for expressions (true if the expression value is >= 0):
    expression ifpositive(expression condexpr, expression trueexpr, expression falseexpr);
    expression andpositive(std::vector<expression> exprs);
    expression orpositive(std::vector<expression> exprs);

    // Evaluate an expression on physical region 'physreg' using interpolation:
    expression on(int physreg, expression expr, bool errorifnotfound = true);
    // Interpolate at coordinates shifted by 'coordshift':
    expression on(int physreg, expression coordshift, expression expr, bool errorifnotfound = true);

    expression comp(int selectedcomp, expression input);
    expression compx(expression input);
    expression compy(expression input);
    expression compz(expression input);

    expression entry(int row, int col, expression input);

    expression transpose(expression input);
    expression inverse(expression input);
    expression determinant(expression input);

    expression grad(expression input);
    expression div(expression input);
    expression curl(expression input);

    // Cross product between vector a and vector b.
    // Any argument that does not have 3 components will be filled with zeros.
    expression crossproduct(expression a, expression b);

    // Frobenius product:
    expression frobeniusproduct(expression a, expression b);

    // Get the trace of a square matrix:
    expression trace(expression a);
    
    // Get the rotation matrix for a 3x3 tensor with (x,y,z) component ordering
    // or a 6x6 tensor in Voigt form. Input angles are in degrees.
    expression rotation(double alphax, double alphay, double alphaz, std::string type = "");

    // Get the determinant of the physical to reference element variable change Jacobian:
    expression detjac(void);

    integration integral(int physreg, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    integration integral(int physreg, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    // For the multiharmonic resolution an extra integer is required
    // to know with how many harmonics the coef multiplying the tf
    // and/or dof should be approximated. Set 'numcoefharms' negative
    // and it will be as if you were calling the above non
    // multiharmonic functions.
    integration integral(int physreg, int numcoefharms, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);
    integration integral(int physreg, int numcoefharms, expression meshdeform, expression tointegrate, int integrationorderdelta = 0, int blocknumber = 0);

    expression dof(expression input, int physreg = -1);
    expression tf(expression input, int physreg = -1);



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
    // Iterative resolution (with or without diagonal scaling):
    void solve(mat A, vec b, vec sol, double& relrestol, int& maxnumit, std::string soltype = "bicgstab", std::string precondtype = "sor", int verbosity = 1, bool diagscaling = false);

    // Generate, solve and save to field a formulation:
    void solve(formulation formul);
    void solve(std::vector<formulation> formuls);
    
    // Convert dB to Nepers:
    expression dbtoneper(expression toconvert);


    ////////// PREDEFINED OPERATORS

    // Gives the engineering strains of a 2D or 3D mechanical displacement vector (Voigt form):
    expression strain(expression input);
    // Gives the Green-Lagrange strains of a 2D or 3D mechanical displacement vector (Voigt form):
    expression greenlagrangestrain(expression gradu);
    // Gives the von Mises stress (Voigt form expected for the argument):
    expression vonmises(expression stress);

    // Weak form of the mass conservation for Navier-Stokes:
    expression predefinedmassconservation(expression dofv, expression tfp, expression rho, expression dtrho, expression gradrho, bool includetimederivs, bool isdensityconstant);
    // Weak form of the inertial forces for Navier-Stokes:
    expression predefinedinertialforce(expression dofv, expression tfv, expression v, expression rho);
    // Weak form of the viscous forces for Navier-Stokes:
    expression predefinedviscousforce(expression dofv, expression tfv, expression mu, bool isdensityconstant, bool isviscosityconstant);


    ////////// PREDEFINED FORMULATIONS

    // Isotropic linear elasticity:
    expression predefinedelasticity(expression dofu, expression tfu, expression Eyoung, expression nupoisson, std::string myoption = "");
    // General anisotropic linear elasticity:
    expression predefinedelasticity(expression dofu, expression tfu, expression elasticitymatrix, std::string myoption = "");

    // Isotropic elasticity with geometrical nonlinearity and prestress (ignored if zero):
    expression predefinedelasticity(expression dofu, expression tfu, field u, expression Eyoung, expression nupoisson, expression prestress, std::string myoption = "");
    // General anisotropic elasticity with geometrical nonlinearity and prestress (ignored if zero):
    expression predefinedelasticity(expression dofu, expression tfu, field u, expression elasticitymatrix, expression prestress, std::string myoption = "");

    expression predefinedelectrostaticforce(expression tfu, expression E, expression epsilon);
    expression predefinedelectrostaticforce(std::vector<expression> dxyztfu, expression E, expression epsilon);
    expression predefinedmagnetostaticforce(expression tfu, expression H, expression mu);
    expression predefinedmagnetostaticforce(std::vector<expression> dxyztfu, expression H, expression mu);

    expression predefinedacoustics(expression dofp, expression tfp, expression soundspeed, expression neperattenuation);
    expression predefinedacousticradiation(expression dofp, expression tfp, expression soundspeed, expression neperattenuation);
    expression predefinedacousticstructureinteraction(expression dofp, expression tfp, expression dofu, expression tfu, expression soundspeed, expression fluiddensity, expression normal, expression neperattenuation, double scaling = 1.0);

    // Stokes flow for Newtonian fluids:
    expression predefinedstokes(expression dofv, expression tfv, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs = false, bool isdensityconstant = true, bool isviscosityconstant = true);
    // Navier-Stokes flow for Newtonian fluids:
    expression predefinednavierstokes(expression dofv, expression tfv, expression v, expression dofp, expression tfp, expression mu, expression rho, expression dtrho, expression gradrho, bool includetimederivs = false, bool isdensityconstant = true, bool isviscosityconstant = true);
};

#endif
