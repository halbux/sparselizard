// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPERATION_H
#define OPERATION_H

#include <iostream>
#include <vector>
#include "densematrix.h"
#include "universe.h"
#include "harmonic.h"
#include "jacobian.h"
#include <cmath>
#include <memory>
#include "rawfield.h"
#include "parameter.h"
#include "selector.h"
#include "elementselector.h"
#include "hierarchicalformfunction.h"
#include "oncontext.h"

// ALL SON-OPERATIONS ARE INCLUDED AT THE END OF THIS HEADER.

class rawfield;
class oncontext;

class operation : public std::enable_shared_from_this<operation>
{

    private:
        
    public:
        
        // Print the operation:
        virtual void print(void) {};
        
        // Interpolate the operation on the evaluation coordinates. 
        // The returned value taken at [harm][0] gives the value of the 
        // harmonic harm. If ...[harm].size() is zero the harmonic is zero.
        // This function can be used for non multiharmonic as well as
        // LINEAR multiharmonic operations. In the first case only 
        // harmonic 1 (cos0) can be used.
        virtual std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        // 'multiharmonicinterpolate' works for general nonlinear 
        // multiharmonic operations. It returns a matrix in which each 
        // column corresponds to a data point (e.g. an operation evaluated 
        // at an evaluation point) and each of the 'numtimeevals' rows to 
        // a time value at which the multiharmonic operation was computed. 
        // 'myfft.h' can be used to get back the harmonics from that matrix.
        virtual densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        // This 'interpolate' is used only in the Jacobian computation.
        virtual std::vector<std::vector<densematrix>> interpolate(int kietaphiderivative, elementselector& elemselect, std::vector<double>& evaluationcoordinates);
        
        virtual bool isdof(void) { return false; };
        virtual bool istf(void) { return false; };
        virtual bool issum(void) { return false; };
        virtual bool isproduct(void) { return false; };
        virtual bool isfield(void) { return false; };
        virtual bool isconstant(void) { return false; };
        
        // True if the expression is a constant 0:
        bool iszero(void);
        
        // True if the expression includes or is a dof/tf:
        virtual bool isdofincluded(void);
        virtual bool istfincluded(void);
        
        // True if the operation only includes harmonic 1:
        virtual bool isharmonicone(std::vector<int> disjregs);
        
        // Get the value of a constant expression:
        virtual double getvalue(void) {};
        
        // Get the field pointer of expressions including a field:
        virtual std::shared_ptr<rawfield> getfieldpointer(void);
        
        // Remove the term of a sum or product:
        virtual void removeterm(int whichterm) {};
        // Count the number of sum or product terms:
        virtual int count(void) {};
        
        // Set a flag on this operation so that when an operation 
        // 'op' including at least once this operation is 
        // interpolated this operation is only computed once 
        // and then reused for all other occurences in 'op'. 
        virtual void reuseit(bool istobereused);
        virtual bool isreused(void);
        
        // Set derivatives for fields, dofs and tfs:
        virtual void setspacederivative(int whichderivative);
        virtual void setkietaphiderivative(int whichderivative);
        virtual void increasetimederivativeorder(int derivativeorder);
        
        // Get info for fields, dofs and tfs:
        virtual int getphysicalregion(void) {};
        virtual int getspacederivative(void) {};
        virtual int gettimederivative(void) {};
        virtual int getkietaphiderivative(void) {};
        
        // Get the 'argnum'th argument:
        virtual std::shared_ptr<operation> getargument(int argnum) {};
        // Replace the 'argnum'th argument:
        virtual void replaceargument(int argnum, std::shared_ptr<operation> newarg) {};
        
        // Get the form function component number used in the field, dof or tf:
        virtual int getformfunctioncomponent(void) {};
        // Know which subfield of the original field it was:
        virtual int getfieldcomponent(void) {};

        // True if the operation can be interpolated on all elements in 
        // the disjoint regions, no matter their total orientation number:
        virtual bool isvalueorientationdependent(std::vector<int> disjregs);
        
        // Duplicate the operation (argument operations are not duplicated!):
        virtual std::shared_ptr<operation> copy(void) {};

        // Get the arguments of the operation (if any):
        virtual std::vector<std::shared_ptr<operation>> getarguments(void) { return {}; };

        // Expand the operation:
        virtual std::shared_ptr<operation> expand(void) { return shared_from_this(); };
        // Group all sum/product terms together:
        virtual void group(void) {};
        // Simplify the operation as it is on the disjoint regions:
        virtual std::shared_ptr<operation> simplify(std::vector<int> disjregs) { return shared_from_this(); };
     
        // Evaluate an operation that only contains x, y and/or z fields without derivatives.
        virtual std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);
        
        // For dof interpolation:
        virtual bool ison(void) {};
        virtual void setoncontext(oncontext& cntxt) {};
        virtual oncontext* getoncontext(void) {};
};

#include "opabs.h"
#include "opacos.h"
#include "opasin.h"
#include "opatan.h"
#include "opathp.h"
#include "opcondition.h"
#include "opconstant.h"
#include "opcos.h"
#include "opcustom.h"
#include "opdetjac.h"
#include "opdof.h"
#include "opestimator.h"
#include "opfield.h"
#include "opfieldorder.h"
#include "opharmonic.h"
#include "opinversion.h"
#include "opinvjac.h"
#include "opjac.h"
#include "oplog10.h"
#include "opmeshsize.h"
#include "opmod.h"
#include "opon.h"
#include "opparameter.h"
#include "oppower.h"
#include "opproduct.h"
#include "opsin.h"
#include "opspline.h"
#include "opsum.h"
#include "optan.h"
#include "optf.h"
#include "optime.h"

#endif
