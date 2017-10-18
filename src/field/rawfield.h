// The 'rawfield' object is the actual field object. The 'field' 
// object is just a user-friendly wrapper for a 'rawfield' object.
//
// A 'rawfield' can include subfields (e.g. for "h1xyz" it includes 3 
// "h1" subfields) in 'mysubfields[i][0]'. It can also include harmonics 
// in 'myharmonics'. In that case 'mysubfields' must be empty. Harmonic
// h is available at myharmonics[h][0] if it exists (otherwise 
// myharmonics[h] is empty).
// When 'mysubfields' AND 'myharmonics' are empty then the raw field
// actually stores field data in its containers.
//
// Note: a raw field can include subfields which can include harmonics but:
//      - an included subfield can not itself include subfields, only harmonics
//      - an included harmonic can not itself include subfields or harmonics

#ifndef RAWFIELD_H
#define RAWFIELD_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "coefmanager.h"
#include "universe.h"
#include "expression.h"
#include "integration.h"
#include "mathop.h"
#include "harmonic.h"
#include "densematrix.h"
#include "element.h"
#include "elements.h"
#include "nodes.h"
#include "disjointregions.h"
#include "hierarchicalformfunction.h"
#include "hierarchicalformfunctioncontainer.h"
#include "hierarchicalformfunctioniterator.h"
#include "lagrangeformfunction.h"
#include "elementselector.h"
#include "memory.h"
#include "selector.h"
#include "field.h"
#include "vectorfieldselect.h"

class vectorfieldselect;
class coefmanager;
class rawvec;
class expression;

class rawfield : public std::enable_shared_from_this<rawfield>
{
    
	private:
        
        bool multiharmonic;

        std::string myname = "";
        std::string mytypename = "";
        
        std::vector<std::vector<   shared_ptr<rawfield>   >> mysubfields = {};
        std::vector<std::vector<   shared_ptr<rawfield>   >> myharmonics = {};
    
        
        // In case there is neither a subfield nor a harmonic, i.e. both 
        // vectors above are empty, then the raw field actually stores 
        // field data in the containers below.
        // In this case the harmonic is said to be equal 1.

		shared_ptr<coefmanager> mycoefmanager = NULL;
        
        // interpolationorder[disjreg] gives the interpolation 
        // order of the field on disjoint region 'disjreg'.
        std::vector<int> interpolationorder = {};
        // myconstraints[disjreg] gives the integration object to compute the 
        // constraint value on the disjoint region. NULL means unconstrained.
        std::vector<shared_ptr<integration>> myconstraints = {};
                
	public:
        
        rawfield(std::string fieldtypename, const std::vector<int> harmonicnumbers, bool ismultiharm);
        rawfield(void) {};
        
        bool ismultiharmonic(void) { return multiharmonic; };
        
        int countcomponents(void);
        int countsubfields(void) { return mysubfields.size(); };
        int countformfunctioncomponents(void);
        
        // List all harmonics in the raw field.
        std::vector<int> getharmonics(void);
        int getfirstharmonic(void);
        bool isharmonicincluded(int harmonic);
        // Print a string showing the harmonics in the field:
        void printharmonics(void);
        
        // Print the raw field name:
        void print(void);
        // Set the raw field name:
        void setname(std::string name) { myname = name; };
        std::string gettypename(void) { return mytypename; };

        void setorder(int physreg, int interpolorder);
        void setconstraint(int physreg, expression input, int extraintegrationdegree = 0);
        // Set homogeneous Dirichlet constraints:
        void setconstraint(int physreg);
        
        shared_ptr<rawfield> getpointer(void) { return shared_from_this(); };

        // Transfer data from a solution vector to the field.
        void getdata(int physreg, vectorfieldselect myvec);
                
        // Select a component.
        shared_ptr<rawfield> comp(int component);
        // Select a single or several harmonics. 
        // Outputs all components corresponding to that harmonic.
        shared_ptr<rawfield> harmonic(int harmonicnumber) { return harmonic(std::vector<int>{harmonicnumber}); };
        shared_ptr<rawfield> harmonic(const std::vector<int> harmonicnumbers);
        
        bool isconstrained(int disjreg) { return not(myconstraints[disjreg] == NULL); };
        std::vector<shared_ptr<integration>> getconstraints(void) { return myconstraints; };
        // Get the interpolation order on a disjoint region.
        // Only valid for fields without subfields.
        int getinterpolationorder(int disjreg);
        
        // Give an error if all harmonics have not the same interpolation order.
        // Only valid for fields without subfields.
        void errornotsameinterpolationorder(int disjreg);
        

        // This interpolate is called in practice:
        std::vector<std::vector<densematrix>> interpolate(int whichderivative, int formfunctioncomponent, elementselector& elemselect, std::vector<double>& evaluationcoordinates);
        
        // The function works only on fields that are not containers.
        densematrix getcoefficients(int elementtypenumber, int interpolorder, std::vector<int> elementnumbers);
        // 'interpolate' outputs the field value at the evaluation coordinates
        // provided as second argument for all elements in 'elementlist'.
        // Set 'whichderivative' to 0, 1, 2 or 3 to get respectively the 
        // no derivative, ki, eta or phi derivative of the field.
        std::vector<std::vector<densematrix>> interpolate(int whichderivative, int formfunctioncomponent, int elementtypenumber, int totalorientation, int interpolorder, std::vector<int> elementnumbers, std::vector<double>& evaluationcoordinates);

};


#endif
