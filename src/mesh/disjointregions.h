// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


#ifndef DISJOINTREGIONS_H
#define DISJOINTREGIONS_H

#include <iostream>
#include <vector>
#include "element.h"

class disjointregions
{

    private:
        
        std::vector<int> rangebegin;
        std::vector<int> rangeend;
        
        // disjointregionsdefinition[i][j] is true if disjoint 
        // region number i is in the jth physical region.
        std::vector<std::vector<bool>> disjointregionsdefinition;
        // 'elementtypenumbers[i]' gives the unique element
        // type number of the elements in the disjoint region i.
        std::vector<int> elementtypenumbers;
        
    public:
        
        int count(void);
        int countelements(int disjointregionnumber);
        
        // Append the disjoint region including elements of a given type number
        // and defined by 'physicalregionsincludingit'. 'physicalregionsincludingit[i]'
        // is true if the ith physical region geometrically covers the disjoint region.
        // The output is the assigned disjoint region number.
        int append(int elementtypenumber, std::vector<bool>& physicalregionsincludingit);
        
        void setrangebegin(int disjointregionnumber, int startrange);
        void setrangeend(int disjointregionnumber, int endrange);
        
        int getrangebegin(int disjointregionnumber);
        int getrangeend(int disjointregionnumber);
        
        int getelementtypenumber(int disjointregionnumber);
        int getelementdimension(int disjointregionnumber);
        
        // Get the numbers of all disjoint regions of a given dimension/element type:
        std::vector<int> getindim(int dim);
        std::vector<int> getintype(int elementtypenumber);
        
        bool isinphysicalregion(int disjointregionnumber, int physicalregionindex);
        
        void removephysicalregions(std::vector<bool> istoremove);
        
        // Clear the content of this object:
        void clear(void);

};

#endif
