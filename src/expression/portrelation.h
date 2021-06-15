// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

// This object stores a port relation. In the multiharmonic case
// the relation for every harmonic is extracted and stored.

#ifndef PORTRELATION_H
#define PORTRELATION_H

#include <memory>
#include "rawport.h"
#include "expression.h"

class portrelation
{

    private:
    
        // Coefficient of each port (cannot be multiharmonic):
        std::vector<expression> mycoefs = {};
        
        // Optional no-port term (cannot be multiharmonic):
        std::vector<expression> mynoportterm = {};
        
    
        // Below are the containers for every relation. Entry [h] gives
        // all terms for the harmonic h relation (no terms if no relation).
        
        // Rawports:
        std::vector< std::vector<std::shared_ptr<rawport>> > myrawports = {};
        
        // Index of each associated coefficient:
        std::vector< std::vector<int> > mycoefinds = {};
    
        // KCM targets:
        std::vector< std::vector<int> > mykcm = {};
    
        // Factors and f0 powers (use 1.0 and 0.0 for none):
        std::vector< std::vector<std::pair<double, int>> > myfactors = {};
    
    public:
    
        portrelation(expression prtrel);
        
        // Count the number of sub-relations:
        int count(void);
        
        // Get all (possibly duplicated) rawports in the relation:
        std::vector<std::shared_ptr<rawport>> getrawports(void);
        
        bool hasnoportterm(void);
        // Evaluate the no-port term
        double evalnoportterm(void);
        
        // For each term in KCM get the associated rawport, relation index and value: 
        void evalrelations(int KCM, std::vector<std::shared_ptr<rawport>>& rps, std::vector<int>& relinds, std::vector<double>& relvals);

};

#endif

