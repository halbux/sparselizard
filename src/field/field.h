// sparselizard - Copyright (C) 2017-2018 A. Halbach and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <alexandre.halbach at ulg.ac.be>.


// This object is a wrapper of the actual field object 'rawfield' pointed
// by 'rawfieldptr' and to which most of the functions are redirected.
// The purpose of this object is to wrap the 'rawfield' object for a 
// convenient user experience.
// An additional advantage is that the shared_ptr type pointer ensures
// the pointed 'rawfield' object is always available when there is
// at least one field using it.
        
#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <vector>
#include <string>
#include "vec.h"
#include "memory.h"
#include "expression.h"
#include "parameter.h"
#include "rawfield.h"
#include "vectorfieldselect.h"

class rawfield;

class field
{
    
	private:
	
        // The actual field:
        shared_ptr<rawfield> rawfieldptr = NULL;
        
	public:

        // Provide the form function type name for the field.
        field(std::string fieldtypename);
        // Also provide the harmonics for a multiharmonic field:
        field(std::string fieldtypename, const std::vector<int> harmonicnumbers);
        // The constructor below should not be used by the user:
        field(shared_ptr<rawfield> rawfieldpointer) { rawfieldptr = rawfieldpointer; };
        
        // Give the number of components in the field.
        int countcomponents(void);
        
        // List all harmonics in the field.
        std::vector<int> getharmonics(void);
        // Print a string showing the harmonics in the field.
        void printharmonics(void);
        
        // Set the field name.
        void setname(std::string name);
        // Print the field name.
        void print(void);
        
        // Set the interpolation order on a physical region.
        void setorder(int physreg, int interpolorder);
        
        // Set a value for the field on a given geometrical region.
        // Use the default order + 'extraintegrationdegree' to 
        // compute the finite element discretisation of 'input'.
        void setvalue(int physreg, expression input, int extraintegrationdegree = 0);
        // Set a zero value:
        void setvalue(int physreg);
        
        // Set an 'input' valued constraint on a physical region. 
        // Use the default order + 'extraintegrationdegree' to 
        // compute the finite element discretisation of 'input'.
        void setconstraint(int physreg, expression input, int extraintegrationdegree = 0);
        // Set an homogeneous Dirichlet constraint.
        void setconstraint(int physreg);
        
        // Transfer data from a field in the solution vector to this field.
        // Only the data corresponding to the physical region is transferred.
        // Transfer data from field 'a' in 'vector|a' to the current field:
        void setdata(int physreg, vectorfieldselect myvec);
        // Transfer data from and to the current field:
        void setdata(int physreg, vec myvec);
        
        shared_ptr<rawfield> getpointer(void) { return rawfieldptr; };
        
        // Select a component.
        field comp(int component);
        field compx(void) { return comp(0); };
        field compy(void) { return comp(1); };
        field compz(void) { return comp(2); };
        
        // Select a single or multiple harmonics. 
        // Outputs all components corresponding to the harmonic.
        field harmonic(int harmonicnumber) { return harmonic(std::vector<int>{harmonicnumber}); };
        field harmonic(const std::vector<int> harmonicnumbers);
        field sin(int freqindex) { return harmonic(2*freqindex); };
        field cos(int freqindex) { return harmonic(2*freqindex+1); };
        
        
        

        
        double integrate(int physreg, expression meshdeform, int integrationorder);
        double integrate(int physreg, int integrationorder);
        
        void write(int physreg, int numfftharms, std::string filename, int lagrangeorder = 1);
        void write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder = 1);

        void write(int physreg, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        void write(int physreg, expression meshdeform, std::string filename, int lagrangeorder = 1, int numtimesteps = -1);
        
        
		// Defining the +, -, * and / operators:
        expression operator+(void);
        expression operator-(void);
        
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
expression operator+(double, field);
expression operator-(double, field);
expression operator*(double, field);
expression operator/(double, field);

expression operator+(parameter&, field);
expression operator-(parameter&, field);
expression operator*(parameter&, field);
expression operator/(parameter&, field);

#endif
