// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <iostream>
#include <string>
#include <vector>
#include "universe.h"
#include "operation.h"
#include "field.h"
#include "disjointregions.h"
#include "physicalregions.h"
#include "physicalregion.h"
#include "elementselector.h"
#include "disjointregionselector.h"
#include <algorithm>
#include "gausspoints.h"
#include "jacobian.h" 
#include <memory>
#include "parameter.h"
#include "polynomial.h"
#include "myfft.h"
#include <cmath>
#include "rawfield.h"
#include "mystring.h"
#include "wallclock.h"
#include "shape.h"
#include "myalgorithm.h"
#include "iointerface.h"
#include "spline.h"
#include "referencecoordinategroup.h"


class vec;
class operation;
class parameter;
class field;
class shape;

class expression
{
    friend class opon;

    private:
        
        int mynumrows = 0;
        int mynumcols = 0;
        
        // myoperations[i*mynumcols+j] gives the 
        // expression operation defined at row i, column j.
        std::vector<std::shared_ptr<operation>> myoperations = {};
        
        std::vector<std::pair<std::string,expression>> inrefcoord = {};
        
        
        // FUNCTIONS TO BE CALLED BY THE PUBLIC FUNCTIONS:
        
        std::vector<double> max(int physreg, expression* meshdeform, int refinement, std::vector<double> xyzrange);
        void interpolate(int physreg, expression* meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        void interpolate(int physreg, expression* meshdeform, std::vector<double>& xyzcoord, std::vector<std::vector<double>>& interpolated, std::vector<bool>& isfound, int numtimeevals);
        double integrate(int physreg, expression* meshdeform, int integrationorder);
        void write(int physreg, int numfftharms, expression* meshdeform, std::string filename, int lagrangeorder, int numtimesteps);
        
    public:
        
        expression(void) {};
        // Implicit conversion from field, double and parameter to expression:
        expression(field);
        expression(double);
        expression(parameter&);
        expression(int numrows, int numcols, std::vector<expression>);
        // Concatenate expressions to create a new one:
        expression(const std::vector<std::vector<expression>> input);
        // Expression whose value depends on if the first argument is greater or equal to zero.
        // If true the expression value is the expression as second argument, if false it is the 
        // expression provided as third argument.
        expression(expression condexpr, expression exprtrue, expression exprfalse);
        // Expression based on a spline interpolation of a discrete function of argument 'arg':
        expression(spline spl, expression arg);
        // Piecewise expression definition:
        expression(std::vector<double> pos, std::vector<expression> exprs, expression tocompare);
        // Custom expression based on a user-defined function:
        expression(int m, int n, std::vector<densematrix> customfct(std::vector<densematrix>), std::vector<expression> exprs);
        
        // Define a 1x1 expression from an operation:
        expression(std::shared_ptr<operation>);
        
        int countrows(void) { return mynumrows; };
        int countcolumns(void) { return mynumcols; };
        
        // Get a given row/column in a matrix expression:
        expression getrow(int rownum);
        expression getcolumn(int colnum);
        
        void reorderrows(std::vector<int> neworder);
        void reordercolumns(std::vector<int> neworder);
        
        // Get the max/min value. All elements will be split 'refinement' times in each direction 
        // to approximate the max/min value and position. Increase 'refinement' for more accuracy.
        // Set {xrangemin, xrangemax, yrangemin, yrangemax, zrangemin, zrangemax} to get the 
        // max/min in a x, y and z bounded domain (optional).
        // The output is {maxvalue, maxxcoord, maxycoord, maxzcoord}.
        std::vector<double> max(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        
        // Interpolate the expression at N (x,y,z) coordinates (provided in 'xyzcoord' in format
        // {x1,y1,z1, x2,y2,z2,...}). After the call the interpolated values of the expression 
        // are in 'interpolated' (non-scalar expressions are flattened and their interpolated 
        // values concatenated one after the other).
        // Only the highest dimension elements in physical region 'physreg' are considered. 
        // In case the ith coordinate is not in the physical region or there was any other 
        // issue then 'isfound[i]' is false. 
        void interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        void interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        // These two functions are added for convenience and work only to interpolate at a single (x,y,z) coordinate.
        // In case the coordinate is not in the physical region or there was any other issue the returned vector is empty.
        std::vector<double> interpolate(int physreg, const std::vector<double> xyzcoord);
        std::vector<double> interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord);
        
        double integrate(int physreg, int integrationorder);
        double integrate(int physreg, expression meshdeform, int integrationorder);
        
        // Compute an FFT transform of the expression and save all the 'numfftharms' harmonics:
        void write(int physreg, int numfftharms, std::string filename, int lagrangeorder = 1);
        void write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder = 1);
        // Save at 'numtimesteps' timesteps. Set -1 to save the harmonics for linear expressions.
        void write(int physreg, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        void write(int physreg, expression meshdeform, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        
        // Save to disk the part of the stream lines that lie on physical region 'physreg'.
        // The stream lines are grown (upstream and downstream) starting from 'startcoords'.
        // 'stepsize' is related to the distance between stream direction updates (decrease for more accuracy).
        void streamline(int physreg, std::string filename, const std::vector<double>& startcoords, double stepsize, bool downstreamonly = false);
        
        // Set a flag on this expression so that when an expression 
        // 'expr' including at least once this expression is 
        // interpolated this expression is only computed once 
        // and then reused for all other occurences in 'expr'. 
        void reuseit(bool istobereused = true);
        
        bool isscalar(void) { return (mynumrows == 1 && mynumcols == 1); };
        bool isharmonicone(std::vector<int> disjregs);
        bool isvalueorientationdependent(std::vector<int> disjregs);
        bool iszero(void);
        
        // Output a vector based on field 'onefield' that stores the barycenter values of the expression.
        vec atbarycenter(int physreg, field onefield);
        
        // Print the expression:
        void print(void);
        
        // Read the function documentation before using it!
        void rotate(double ax, double ay, double az, std::string leftop = "default", std::string rightop = "default");
        
        expression at(int row, int col);
        
        // Evaluate a scalar expression that only contains x, y and/or z fields without derivatives.
        std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);
        
        // Output the resized expression (filled with zero if larger):
        expression resize(int numrows, int numcols);
        
        
        
        // THE FUNCTIONS BELOW ARE NOT MEANT TO BE CALLED BY THE USER!
        
        // 'whichderivative' is 1 for x, 2 for y and 3 for z.
        expression spacederivative(int whichderivative);
        // Now the derivative in the reference element.
        // 'whichderivative' is 1 for ki, 2 for eta and 3 for phi.
        expression kietaphiderivative(int whichderivative);
        expression timederivative(int derivativeorder);
        
        std::vector<std::pair<std::string,expression>> getinrefcoord(void);
        
        expression transpose(void);
        // Get the submatrix obtained by removing a row and a column:
        expression removerowandcol(int rowtoremove, int coltoremove);
        expression determinant(void);
        expression cofactormatrix(void);
        expression invert(void);
        
        expression pow(expression);
        // In dof() and tf() '-1' selects no physical region:
        expression dof(int physreg);
        expression tf(int physreg);
        expression sin(void);
        expression cos(void);
        expression tan(void);
        expression asin(void);
        expression acos(void);
        expression atan(void);
        expression abs(void);
        expression log10(void);
        expression mod(double modval);
        
        expression on(int physreg, expression* coordshift, bool errorifnotfound);
        
        // The time variable:
        expression time(void);
        
        expression invjac(int row, int col);
        expression jac(int row, int col);
        expression detjac(void);
        // Get the whole 3x3 Jacobian matrix:
        expression invjac(void);
        expression jac(void);
        
        expression getcopy(void);
        
        std::shared_ptr<operation> getoperationinarray(int row, int col);
        
        // Expand only the terms containing a dof or a testfun.
        // The expression must be scalar!
        void expand(void);
        // After having expanded the expression its form should be a sum
        // of terms like coef*dof*tf where coef does not include a dof or 
        // tf and where different operations can be applied to the dof and 
        // tf from one term in the sum to another. All coefs are put in
        // the output at output[0], all dofs in output[1] and all tfs 
        // in output[2]. A coef can then be found at output[0][s][i]. 
        // All coefs corresponding to dof-tf pairs with the same applied 
        // space and time derivatives are added together in a new operation.
        // They are added together at an a priori unknown index i. 
        // output[0], output[1] and output[2] may have multiple slices 
        // output[0][s] each corresponding to a unique dof field*-tf field* 
        // pair (with a unique combination of applied time derivatives). 
        //
        // The expression must be scalar!
        std::vector< std::vector<std::vector<std::shared_ptr<operation>>> > extractdoftfpolynomial(int elementdimension);
        
        
        // Defining the +, -, * and / operators:
        expression operator+(void);
        expression operator-(void);
        
        expression operator+(expression);
        expression operator-(expression);
        expression operator*(expression);
        expression operator/(expression);
        
        expression operator+(field);
        expression operator-(field);
        expression operator*(field);
        expression operator/(field);
        
        expression operator+(double);
        expression operator-(double);
        expression operator*(double);
        expression operator/(double);
        
        expression operator+(parameter&);
        expression operator-(parameter&);
        expression operator*(parameter&);
        expression operator/(parameter&);
};

// Define the left version of the operators based on the right one.
expression operator+(double, expression);
expression operator-(double, expression);
expression operator*(double, expression);
expression operator/(double, expression);

expression operator+(field, expression);
expression operator-(field, expression);
expression operator*(field, expression);
expression operator/(field, expression);

expression operator+(parameter&, expression);
expression operator-(parameter&, expression);
expression operator*(parameter&, expression);
expression operator/(parameter&, expression);

#endif



