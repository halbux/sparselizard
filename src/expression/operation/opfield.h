// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef OPFIELD_H
#define OPFIELD_H

#include "operation.h"

class opfield: public operation
{

    private:

        // Fields are always reused

        int timederivativeorder = 0;
        // 0 is no derivative, 1 is x, 2 is y and 3 is z.
        int spacederivative = 0;
        // 0 is no derivative, 1 is ki, 2 is eta and 3 is phi.
        int kietaphiderivative = 0;

        // The selected form function component in the field:
        int formfunctioncomponent = 0;
        // For printing purposes:
        int fieldcomponent = -1;

        std::shared_ptr<rawfield> myfield;

    public:

        opfield(std::shared_ptr<rawfield> fieldin) { myfield = fieldin; };

        void setspacederivative(int whichderivative);
        void setkietaphiderivative(int whichderivative);
        void increasetimederivativeorder(int amount);
        void selectformfunctioncomponent(int comp) { formfunctioncomponent = comp; };

        void setfieldcomponent(int comp) { fieldcomponent = comp; };

        bool isfield(void) { return true; };

        bool isharmonicone(std::vector<int> disjregs);

        // This 'interpolate' is used only in the Jacobian computation.
        std::vector<std::vector<densematrix>> interpolate(int kietaphiderivative, elementselector& elemselect, std::vector<double>& evaluationcoordinates);

        std::vector<std::vector<densematrix>> interpolate(elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);
        densematrix multiharmonicinterpolate(int numtimeevals, elementselector& elemselect, std::vector<double>& evaluationcoordinates, expression* meshdeform);

        std::shared_ptr<rawfield> getfieldpointer(void) { return myfield; };
        int getformfunctioncomponent(void) { return formfunctioncomponent; };
        int getfieldcomponent(void) { return fieldcomponent; };
        int getspacederivative(void) { return spacederivative; };
        int getkietaphiderivative(void) { return kietaphiderivative; };
        int gettimederivative(void) { return timederivativeorder; };

        bool isvalueorientationdependent(std::vector<int> disjregs);

        std::shared_ptr<operation> copy(void);

        std::vector<double> evaluate(std::vector<double>& xcoords, std::vector<double>& ycoords, std::vector<double>& zcoords);

        void print(void);

};

#endif
