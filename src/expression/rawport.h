// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// A rawport can either be associated to a rawfield or live its own life.

#ifndef RAWPORT_H
#define RAWPORT_H

#include "rawfield.h"

class rawfield;

class rawport : public std::enable_shared_from_this<rawport>
{

    private:
    
        bool amimultiharmonic = false;
        
        std::string myname = "";
        
        std::vector<std::vector<   std::shared_ptr<rawport>   >> myharmonics = {};
        
        // The data below is not defined if the above vector is not empty.
        
        bool myisprimal = true;
        
        // Dual/primal port if this is the primal/dual.
        std::shared_ptr<rawport> mybrother = NULL;
        
        // This port is associated to a unique rawfield and physical region:
        int myphysreg = -1;
        std::shared_ptr<rawfield> myrawfield = NULL;
    
    public:
    
        rawport(void) {};
        rawport(std::vector<int> harmonicnumbers, bool ismultiharm);
    
        bool ismultiharmonic(void) { return amimultiharmonic; };
    
        void setname(std::string name);
        std::string getname(void);
        
        std::vector<int> getharmonics(void);
        std::shared_ptr<rawport> harmonic(int harmonicnumber);
        std::shared_ptr<rawport> harmonic(std::vector<int> harmonicnumbers);
        
        bool isprimal(void);
        
        int getphysicalregion(void);
        
        std::shared_ptr<rawfield> getrawfield(void);
        
        std::shared_ptr<rawport> getbrother(void);
        
        // Check if this rawport is associated to a rawfield:
        bool isassociated(void);
        
        // Associate the rawport to a rawfield:
        void associate(bool isprim, std::shared_ptr<rawport> bro, int physreg, std::shared_ptr<rawfield> rf);
        
};

#endif

