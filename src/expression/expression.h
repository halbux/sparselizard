// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <iostream>
#include <string>
#include <vector>
#include "universe.h"
#include "operation.h"
#include "field.h"
#include "disjointregions.h"
#include "gmshinterface.h"
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


class vec;
class operation;
class parameter;
class field;

class expression
{    
	private:
                
        int mynumrows = 1;
        int mynumcols = 1;
        
        // myoperations[i*mynumcols+j] gives the 
        // expression operation defined at row i, column j.
        std::vector<std::shared_ptr<operation>> myoperations;
        
        // Is this expression the projection on the physical element of a 
        // vector form function based-field defined in the reference element?
        // If no then 'unprojectedfield' has size zero.
        std::vector<expression> unprojectedfield = {};
        
        
        // FUNCTIONS TO BE CALLED BY THE PUBLIC FUNCTIONS:
        
        std::vector<double> max(int physreg, expression* meshdeform, int refinement, std::vector<double> xyzrange);
        double integrate(int physreg, expression* meshdeform, int integrationorder);
        void write(int physreg, int numfftharms, expression* meshdeform, std::string filename, int lagrangeorder, int numtimesteps);
        
        
	public:
        
        expression(void) {};
        // Implicit conversion from field, double and parameter to expression:
        expression(field);
        expression(double);
        expression(parameter&);
        expression(int numrows, int numcols, std::vector<expression>);
        
        int countrows(void) { return mynumrows; };
        int countcolumns(void) { return mynumcols; };

        // Get the max/min value. All elements will be split 'refinement' times in each direction 
        // to approximate the max/min value and position. Increase 'refinement' for more accuracy.
		// Set {xrangemin, xrangemax, yrangemin, yrangemax, zrangemin, zrangemax} to get the 
		// max/min in a x, y and/or z bounded domain (optional).
        // The output is {maxvalue, maxxcoord, maxycoord, maxzcoord}.
        std::vector<double> max(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});

        double integrate(int physreg, int integrationorder);
        double integrate(int physreg, expression meshdeform, int integrationorder);
        
        // Compute an FFT transform of the expression and save all the 'numfftharms' harmonics:
        void write(int physreg, int numfftharms, std::string filename, int lagrangeorder = 1);
        void write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder = 1);
        // Save at 'numtimesteps' timesteps. Set -1 to save the harmonics for linear expressions.
        void write(int physreg, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        void write(int physreg, expression meshdeform, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        
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
        // Output a vector based on field 'onefield' that stores the value of the 
        // expression integrated on every element separately in region 'physreg'.
        vec integrateonelements(int physreg, field onefield, int integrationorder);
        
        // Print the expression:
        void print(void);
        
        expression at(int row, int col);
           
        
        // THE FUNCTIONS BELOW ARE NOT MEANT TO BE CALLED BY THE USER!
        
        // 'whichderivative' is 1 for x, 2 for y and 3 for z.
        expression spacederivative(int whichderivative);
        // Now the derivative in the reference element.
        // 'whichderivative' is 1 for ki, 2 for eta and 3 for phi.
        expression kietaphiderivative(int whichderivative);
        expression timederivative(int derivativeorder);
        
        bool isprojectedfield(void) { return (unprojectedfield.size() > 0); };
        expression getunprojectedfield(void);
        
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
        expression abs(void);
        expression log10(void);
        
        // The time variable:
        expression time(void);
        
        expression invjac(int row, int col);
        expression jac(int row, int col);
        expression detjac(void);
        // Get the whole 3x3 matrix (with zeros everywhere but on the 
        // 1x1 submatrix in 1D and on the 2x2 submatrix in 2D). 
        expression invjac(void);
        expression jac(void);
        
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



