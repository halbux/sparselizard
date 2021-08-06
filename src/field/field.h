// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object is a wrapper of the actual field object 'rawfield' pointed
// by 'rawfieldptr' and to which most of the functions are redirected.
// The purpose of this object is to wrap the 'rawfield' object for a 
// convenient user experience.
// An additional advantage is that the std::shared_ptr type pointer ensures
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
#include "spanningtree.h"
#include "port.h"

class spanningtree;
class vectorfieldselect;
class parameter;
class rawfield;
class port;

class field
{

    private:

        // The actual field:
        std::shared_ptr<rawfield> rawfieldptr = NULL;
        
        void errorifpointerisnull(void);

    public:

        field(void) {};
        // Provide the form function type name for the field.
        field(std::string fieldtypename);
        // Also provide the harmonics for a multiharmonic field:
        field(std::string fieldtypename, const std::vector<int> harmonicnumbers);
        // Also provide the spanning tree used for gauging fields:
        field(std::string fieldtypename, spanningtree spantree);
        field(std::string fieldtypename, const std::vector<int> harmonicnumbers, spanningtree spantree);
        // The constructor below should not be used by the user:
        field(std::shared_ptr<rawfield> rawfieldpointer) { rawfieldptr = rawfieldpointer; };

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
        void setorder(expression criterion, int loworder, int highorder);
        void setorder(double targeterror, int loworder, int highorder, double absthres);
        
        // Associate the primal and dual port to the field.
        void setport(int physreg, port primal, port dual);

        // Set a value for the field on a given geometrical region.
        // Use the default order + 'extraintegrationdegree' to 
        // compute the finite element discretisation of 'input'.
        void setvalue(int physreg, expression input, int extraintegrationdegree = 0);
        // The 'input' expression is evaluated on the mesh deformed by 'meshdeform':
        void setvalue(int physreg, expression meshdeform, expression input, int extraintegrationdegree = 0);
        // An FFT is used to project the 'input' expression:
        void setvalue(int physreg, int numfftharms, expression input, int extraintegrationdegree = 0);
        void setvalue(int physreg, int numfftharms, expression meshdeform, expression input, int extraintegrationdegree = 0);
        // Set a zero value:
        void setvalue(int physreg);
        
        // Set/get value at nodes for 'h1' type fields:
        void setnodalvalues(intdensematrix nodenumbers, densematrix values);
        densematrix getnodalvalues(intdensematrix nodenumbers);

        // Set an 'input' valued constraint on a physical region. 
        // Use the default order + 'extraintegrationdegree' to 
        // compute the finite element discretisation of 'input'.
        void setconstraint(int physreg, expression input, int extraintegrationdegree = 0);
        // The 'input' expression is evaluated on the mesh deformed by 'meshdeform':
        void setconstraint(int physreg, expression meshdeform, expression input, int extraintegrationdegree = 0);
        // An FFT is used to project the 'input' expression:
        void setconstraint(int physreg, int numfftharms, expression input, int extraintegrationdegree = 0);
        void setconstraint(int physreg, int numfftharms, expression meshdeform, expression input, int extraintegrationdegree = 0);
        // Set an homogeneous Dirichlet constraint.
        void setconstraint(int physreg);

        // Set a 'valexpr' valued constraint on the NODE-ASSOCIATED degrees of freedom of 
        // physical region 'physreg' for which 'condexpr' is greater or equal to zero.
        // All dofs that are not associated to nodes or the nodes at which condexpr < 0 are 
        // left unconstrained (unless constrained by another constraint).
        //
        // This function should only be used for fields with nodal shape functions ("h1" like).
        void setconditionalconstraint(int physreg, expression condexpr, expression valexpr);

        // Set a gauge condition on a given physical region:
        void setgauge(int physreg);

        // Transfer data from a field in the solution vector to this field.
        // Only the data corresponding to the physical region is transferred.
        // 'op' can be 'add' or 'set'.
        // Transfer data from field 'a' in 'vector|a' to the current field:
        void setdata(int physreg, vectorfieldselect myvec, std::string op = "set");
        // Transfer data from and to the current field:
        void setdata(int physreg, vec myvec, std::string op = "set");
        
        // Set the source value at every cut:
        void setcohomologysources(std::vector<int> cutphysregs, std::vector<double> cutvalues);

        // Allow/forbid automatic updating of the field value during hp-adaptivity:
        void automaticupdate(bool updateit);
        void noautomaticupdate(void);
        
        void setupdateaccuracy(int extraintegrationorder);
        
        std::shared_ptr<rawfield> getpointer(void) { return rawfieldptr; };

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



        vec atbarycenter(int physreg, field onefield);

        std::vector<double> max(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> max(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, int refinement, std::vector<double> xyzrange = {});
        std::vector<double> min(int physreg, expression meshdeform, int refinement, std::vector<double> xyzrange = {});
        
        void interpolate(int physreg, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        void interpolate(int physreg, expression meshdeform, std::vector<double>& xyzcoord, std::vector<double>& interpolated, std::vector<bool>& isfound);
        
        std::vector<double> interpolate(int physreg, const std::vector<double> xyzcoord);
        std::vector<double> interpolate(int physreg, expression meshdeform, const std::vector<double> xyzcoord);
        
        double integrate(int physreg, expression meshdeform, int integrationorder);
        double integrate(int physreg, int integrationorder);

        void write(int physreg, int numfftharms, std::string filename, int lagrangeorder);
        void write(int physreg, int numfftharms, expression meshdeform, std::string filename, int lagrangeorder);

        void write(int physreg, std::string filename, int lagrangeorder, int numtimesteps = -1);
        void write(int physreg, expression meshdeform, std::string filename, int lagrangeorder, int numtimesteps = -1);
        
        // Write/load the raw field data to/from compact sparselizard format:
        void writeraw(int physreg, std::string filename, bool isbinary = false, std::vector<double> extradata = {});
        std::vector<double> loadraw(std::string filename, bool isbinary = false);
        

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

        expression operator+(parameter);
        expression operator-(parameter);
        expression operator*(parameter);
        expression operator/(parameter);

};

// Define the left version of the operators based on the right one.
expression operator+(double, field);
expression operator-(double, field);
expression operator*(double, field);
expression operator/(double, field);

expression operator+(parameter, field);
expression operator-(parameter, field);
expression operator*(parameter, field);
expression operator/(parameter, field);

#endif

