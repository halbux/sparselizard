// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// Hierarchical form functions are ordered as follows:
//
// 1. vertex based
// 2. edge based
// 3. face based
// 4. volume based
//
// In each category the form functions are then sorted according to the 
// node, edge, face or volume on which they are defined. E.g. for the 
// vertex based form functions we get first the form functions associated 
// to node 1 then node 2, node 3,...
// On a given node, edge, face or volume the form functions are finally
// sorted with increasing order.

#ifndef HIERARCHICALFORMFUNCTIONITERATOR_H
#define HIERARCHICALFORMFUNCTIONITERATOR_H

#include <iostream>
#include <string>
#include <vector>
#include "element.h"
#include "hierarchicalformfunction.h"
#include <memory>
#include "selector.h"


class hierarchicalformfunctioniterator
{

    private:
    
        std::shared_ptr<hierarchicalformfunction> myformfunction;
        
        int myelementtypenumber;
        int myorder;
        
        // currentdimension is 0 for vertex based, 1 for edge, 2 for face and 3 for volume:
        int currentdimension = 0;
        // Which node/line/face/volume are we at?
        int currentnodeedgefacevolumeindex = 0;
        // The current form function index in the current node/edge/face/volume:
        int currentformfunctionindexinnodeedgefacevolume = -1;
        // The overall form function number:
        int overallformfunctionindex = -1;
        
    public:

        hierarchicalformfunctioniterator(std::string formfunctiontypename, int elementtypenumber, int order);

        // Get the total number of form functions:
        int count(void);
        // 'next' moves to the next defined form function. 
        void next(void);

        int getdimension(void) { return currentdimension; };
        int getnodeedgefacevolumeindex(void) { return currentnodeedgefacevolumeindex; };
        int getformfunctionindexinnodeedgefacevolume(void) { return currentformfunctionindexinnodeedgefacevolume; };
        int getformfunctionindexincurrentorderinnodeedgefacevolume(void);
        int getformfunctionorder(void);
        int getassociatedelementtype(void);
        int getoverallformfunctionindex(void) { return overallformfunctionindex; };
        
        // Print info for the current form function.
        void print(void);
        
};

#endif
