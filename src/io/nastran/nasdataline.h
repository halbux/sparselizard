// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.

#ifndef NASDATALINE_H
#define NASDATALINE_H

#include <iostream>
#include <vector>
#include <string>
#include "logs.h"
#include "mystring.h"

class nasdataline
{
    private:

        // Format 1 for 'small', 2 for 'large' and 3 for 'free'.
        int lineformat = 0;

        bool isitgriddata = false;
        bool isitelementdata = false;

        // For grid data:        
        std::vector<double> nodecoordinate = std::vector<double>(3,0);

        // For element data:
        int elementtypenumber = -1;
        int groupnumber = -1;
        int currentvertexindex = 0;
        std::vector<int> vertices = {};
        
        // Convert the element name to our type number and number of vertices:
        std::vector<int> translateelementname(std::string elemname);
    
    public:

        // Add a data line. True is returned in case any data can be returned
        // (i.e. when at the end of multiple lines, not in a comment,...).
        bool addline(std::string linetoadd);    

        // Does the data correspond to grid data:
        bool isgriddata(void) { return isitgriddata; };
        // Does the data correspond to element data (including the grid number list):
        bool iselementdata(void) { return isitelementdata; };

        ///// FOR GRID DATA. 
        // Only works if data is node data and when 'addline' has returned true.
        std::vector<double> getnodecoordinates(void) { return nodecoordinate; };

        ///// FOR ELEMENT DATA. 
        // Only works if data is element data and when 'addline' has returned true.
        //
        // Get the element type number as defined in sparselizard:
        int getelementtypenumber(void) { return elementtypenumber; };
        int getgroupnumber(void) { return groupnumber; };
        std::vector<int> getvertices(void) { return vertices; };

        // Format 1 for 'small', 2 for 'large' and 3 for 'free'.
        int getformat(void) { return lineformat; };
};

#endif
 
