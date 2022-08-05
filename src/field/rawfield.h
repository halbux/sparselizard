// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


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
#include "sl.h"
#include "harmonic.h"
#include "densemat.h"
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
#include "rawspanningtree.h"
#include "rawmesh.h"
#include "rawport.h"

class rawmesh;
class vectorfieldselect;
class coefmanager;
class rawvec;
class expression;
class elementselector;
class rawspanningtree;

class rawfield : public std::enable_shared_from_this<rawfield>
{
    
    private:
        
        bool amimultiharmonic = false;

        std::string myname = "";
        std::string mytypename = "";
        
        std::vector<std::vector<   std::shared_ptr<rawfield>   >> mysubfields = {};
        std::vector<std::vector<   std::shared_ptr<rawfield>   >> myharmonics = {};
    
        
        // In case there is neither a subfield nor a harmonic, i.e. both 
        // vectors above are empty, then the raw field actually stores 
        // field data in the containers below.
        // In this case the harmonic is said to be equal 1.

        std::shared_ptr<coefmanager> mycoefmanager = NULL;
        
        // interpolationorder[disjreg] gives the interpolation 
        // order of the field on disjoint region 'disjreg'.
        std::vector<int> interpolationorder = {};
        // mydisjregconstraints[disjreg] gives the information to compute the 
        // constraint value on the disjoint region. NULL means unconstrained.
        std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>> mydisjregconstraints = {};
        // myconditionalconstraints[disjreg] gives the {conditional expression, constraint value} 
        // on the NODAL disjoint region 'disjreg'. Empty means unconstrained.
        std::vector<std::vector<expression>> myconditionalconstraints = {};
    
        // The spanning tree used for gauging fields (NULL if none):
        std::shared_ptr<rawspanningtree> myspanningtree = NULL;

        // isitgauged[disjreg] is true if disjoint region 'disjreg' is gauged.
        std::vector<bool> isitgauged = {};
        
        // isitported[disjreg] is true if a port is associated to the disjoint region.
        std::vector<bool> isitported = {};
        
        
        
        bool ispadaptive = false;
        
        std::shared_ptr<ptracker> myptracker = NULL;

        // Track the calls to 'setorder', 'setdisjregconstraint', 'setconditionalconstraint', 'setgauge', 'setport'.
        std::vector<std::pair<int, int>> myordertracker = {};
        std::vector<std::tuple<int, int, std::vector<expression>, expression, int>> mydisjregconstrainttracker = {};
        std::vector<std::tuple<int, expression, expression>> myconditionalconstrainttracker = {};
        std::vector<int> mygaugetracker = {};
        std::vector<std::tuple<int, std::shared_ptr<rawport>, std::shared_ptr<rawport>>> myporttracker = {};
        
        // To avoid infinite recursive calls:
        bool issynchronizing = false;
        // Allow/forbid syncing:
        bool issynchronizingallowed = true;
        // Allow/forbid value syncing:
        bool isvaluesynchronizingallowed = true;
        
        int myupdateaccuracy = 0;
        
        
        // Mesh on which this object is based:
        std::shared_ptr<rawmesh> myrawmesh = NULL;

    public:
  
        // Synchronize with the hp-adapted mesh.
        // If provided, 'physregsfororder' must be {physreg1,orderpr1,physreg2,orderpr2,...} with ORDERS SORTED ASCENDINGLY.
        void synchronize(std::vector<int> physregsfororder = {}, std::vector<int> disjregsfororder = {});
        
        void updateshapefunctions(expression updateexpr, expression* meshdeform, std::vector<std::vector<int>> drsindims, int updateaccuracy, bool withtiming = false);
        void updatenodalshapefunctions(expression updateexpr, expression* meshdeform, std::vector<std::vector<int>> drsindims);
        void updateothershapefunctions(int dim, expression updateexpr, expression* meshdeform, std::vector<std::vector<int>> drsindims, int updateaccuracy);
        
        void allowsynchronizing(bool allowit);
        void allowvaluesynchronizing(bool allowit);
        
        void setupdateaccuracy(int extraintegrationorder);
        
        rawfield(std::string fieldtypename, const std::vector<int> harmonicnumbers, bool ismultiharm);
        rawfield(void);
        // Get a new field with interpolation orders from 'dm' (the corresponding field must have been selected in 'dm'):
        rawfield(dofmanager* dm, std::shared_ptr<rawmesh> rm, std::shared_ptr<ptracker> pt);
        
        ~rawfield(void);
        
        bool ismultiharmonic(void) { return amimultiharmonic; };
        
        int countcomponents(void);
        int countsubfields(void) { return mysubfields.size(); };
        int countformfunctioncomponents(void);
        
        // List all harmonics in the raw field.
        std::vector<int> getharmonics(void);
        int getfirstharmonic(void);
        bool isharmonicincluded(int harmonic);
        // Print a string showing the harmonics in the field:
        void printharmonics(void);
        
        // Reset the coef manager to all zero and return the current one:
        std::shared_ptr<coefmanager> resetcoefmanager(void);
        void setcoefmanager(std::shared_ptr<coefmanager> cm);
        
        // Print the raw field name:
        void print(void);
        void printvalues(bool databoundsonly = true);
        // Set the raw field name:
        void setname(std::string name);
        std::string gettypename(bool familyonly = true);

        void setorder(int physreg, int interpolorder, bool iscalledbyuser = true);
        void setorder(expression criterion, int loworder, int highorder, double critrange); // critrange -1 for automatic choice
        
        // Associate the primal and dual port to the field:
        void setport(int physreg, std::shared_ptr<rawport> primal, std::shared_ptr<rawport> dual);
        
        void setvalue(int physreg, int numfftharms, expression* meshdeform, expression input, int extraintegrationdegree = 0);
        // Set a zero value:
        void setvalue(int physreg);
        
        // Selected elements can include multiple orientations and field orders.
        // 'gpcoordsin' must correspond to Gauss coordinates for the selected element type.
        // If reuse is allowed then all arguments must correspond to the reused data.
        void setvalue(elementselector& elemselect, std::vector<double>& gpcoordsin, expression* meshdeform, densemat values);
        
        // Set/get value at nodes for 'h1' type fields:
        void setnodalvalues(indexmat nodenumbers, densemat values);
        densemat getnodalvalues(indexmat nodenumbers);
        
        // Zero all dofs on the requested physical region:
        void setzerovalue(int physreg);
        
        void setdisjregconstraint(int physreg, int numfftharms, expression* meshdeform, expression input, int extraintegrationdegree = 0);
        // Set homogeneous Dirichlet constraints:
        void setdisjregconstraint(int physreg);
        
        // Set a conditional constraint:
        void setconditionalconstraint(int physreg, expression condexpr, expression valexpr);

        // Set a gauge condition:
        void setgauge(int physreg);

        void setspanningtree(std::shared_ptr<rawspanningtree> spantree);
        // This should only be called on a field without subfields or harmonics:
        std::shared_ptr<rawspanningtree> getspanningtree(void);
        
        std::shared_ptr<rawfield> getpointer(void);
        std::shared_ptr<rawmesh> getrawmesh(void);
        std::shared_ptr<ptracker> getptracker(void);
        
        std::shared_ptr<coefmanager> getcoefmanager(void);

        // Transfer data from a solution vector to the field.
        // Get from all regions with physreg set to -1. 'op' can be 'add' or 'set'. 
        void setdata(int physreg, vectorfieldselect myvec, std::string op = "set");
        
        // Transfer data from the rawfield to a vectorfieldselect:
        void transferdata(int physreg, vectorfieldselect myvec, std::string op);
        
        // Set the source value at every cut:
        void setcohomologysources(std::vector<int> cutphysregs, std::vector<double> cutvalues);
        
        // Select a component.
        std::shared_ptr<rawfield> comp(int component);
        // Select a single or several harmonics. 
        // Outputs all components corresponding to that harmonic.
        std::shared_ptr<rawfield> harmonic(int harmonicnumber);
        std::shared_ptr<rawfield> harmonic(const std::vector<int> harmonicnumbers);
        
        // Get all included rawfields (subfields and harmonics in each subfield):
        std::vector<std::shared_ptr<rawfield>> getsons(void);
        
        // Only valid for fields without subfields.
        bool isdisjregconstrained(int disjreg);
        std::vector<std::shared_ptr<std::tuple<int, int, std::vector<expression>, expression, int, int>>> getdisjregconstraints(void);
        
        bool isconditionallyconstrained(int disjreg);
        std::vector<std::vector<expression>> getconditionalconstraints(void);

        bool isgauged(int disjreg);

        // Get the interpolation order on a disjoint region.
        // Only valid for fields without subfields.
        int getinterpolationorder(int disjreg);
        std::vector<int> getinterpolationorders(void);
        // This function returns the field interpolation order on all elements requested and the max order encountered.
        int getinterpolationorders(int elementtypenumber, std::vector<int>& elementnumbers, std::vector<int>& fieldorders);
        // This function returns the lowest order containing alpha % of the shape function coefficient weight.
        void getinterpolationorders(int fieldorder, double alpha, double absthres, std::vector<double>& weightsforeachorder, std::vector<int>& lowestorders);
        
        // 'weightsforeachorder' has size numelems x fieldorder+1:
        void getweightsforeachorder(int elementtypenumber, int fieldorder, std::vector<int>& elementnumbers, std::vector<double>& weightsforeachorder);
        
        // Give an error if all harmonics have not the same interpolation order.
        // Only valid for fields without subfields.
        void errornotsameinterpolationorder(int disjreg);
        
        // Get the average of the coefficients for all shape functions up to a given order:
        void getaverage(int elementtypenumber, std::vector<int>& elementnumbers, int maxorder, std::vector<double>& averagevals);
        

        // Get a vector listing all subfields and every harmonic for every subfield 
        // as a pair of {subfieldnum,harmnum} and the rawfield pointer:
        std::vector<std::pair<std::vector<int>, std::shared_ptr<rawfield>>> getallsons(void);

        // Write/load the raw data to/from compact sparselizard format:
        void writeraw(int physreg, std::string filename, bool isbinary, std::vector<double> extradata);
        std::vector<double> loadraw(std::string filename, bool isbinary);
        

        // Return {dkix,dkiy,...,detax,detay,...}:
        std::vector<densemat> getjacterms(elementselector& elemselect, std::vector<double>& evaluationcoordinates);


        // This interpolate is called in practice:
        std::vector<std::vector<densemat>> interpolate(int whichderivative, int formfunctioncomponent, elementselector& elemselect, std::vector<double>& evaluationcoordinates);
        
        // The function works only on fields that are not containers.
        densemat getcoefficients(int elementtypenumber, int interpolorder, std::vector<int> elementnumbers);
        // 'interpolate' outputs the field value at the evaluation coordinates
        // provided as second argument for all elements in 'elementlist'.
        // Set 'whichderivative' to 0, 1, 2 or 3 to get respectively the 
        // no derivative, ki, eta or phi derivative of the field.
        std::vector<std::vector<densemat>> interpolate(int whichderivative, int formfunctioncomponent, int elementtypenumber, int totalorientation, int interpolorder, std::vector<int> elementnumbers, std::vector<double>& evaluationcoordinates);

};


#endif
