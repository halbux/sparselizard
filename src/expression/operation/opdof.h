// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPDOF_H
#define OPDOF_H

#include "operation.h"

class opdof: public operation
{

    private:
        
        int timederivativeorder = 0;
        // The space derivative is only temporary. It is split later on into 
        // a product of invjac terms and ki, eta and phi field derivatives.
        // 0 is no derivative, 1 is x, 2 is y and 3 is z.
        int spacederivative = 0;
        // 0 is no derivative, 1 is ki, 2 is eta and 3 is phi.
        int kietaphiderivative = 0;
        
        // The selected form function component in the field:
        int formfunctioncomponent = 0;
        // For printing purposes:
        int fieldcomponent = -1;
        
        std::shared_ptr<rawfield> myfield;
        
        int myphysicalregion;
        
        
        // In case of dof interpolation:
        std::vector<oncontext> myoncontext = {};
    
    public:
        
        opdof(std::shared_ptr<rawfield> fieldin, int physreg = -1);
        
        void setspacederivative(int whichderivative);
        void increasetimederivativeorder(int amount);
        void selectformfunctioncomponent(int comp) { formfunctioncomponent = comp; };

        void setfieldcomponent(int comp) { fieldcomponent = comp; };
        
        void setkietaphiderivative(int whichderivative) { kietaphiderivative = whichderivative; spacederivative = 0; };
        int getkietaphiderivative(void) { return kietaphiderivative; };
        
        bool isdof(void) { return true; };
        
        bool isharmonicone(std::vector<int> disjregs);
        
        std::shared_ptr<rawfield> getfieldpointer(void) { return myfield; };
        bool isphysicalregiondefined(void) { return (myphysicalregion != -1); };
        int getphysicalregion(void) { return myphysicalregion; };
        int getspacederivative(void) { return spacederivative; };
        int gettimederivative(void) { return timederivativeorder; };
        int getformfunctioncomponent(void) { return formfunctioncomponent; };
        int getfieldcomponent(void) { return fieldcomponent; };
        
        bool isvalueorientationdependent(std::vector<int> disjregs);

        std::shared_ptr<operation> copy(void);
        
        void print(void);
        
        
        // In case of dof interpolation:
        bool ison(void);
        void setoncontext(oncontext& cntxt);
        oncontext* getoncontext(void);

};

#endif
