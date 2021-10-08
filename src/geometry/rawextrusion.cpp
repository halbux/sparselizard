#include "rawextrusion.h"


rawextrusion::rawextrusion(int physreg, std::shared_ptr<rawshape> innerrawshape, double height, int numlayers, std::vector<double> extrudedirection)
{
    myphysicalregion = physreg;

    // Sons will be created while meshing

    myheight = height;

    mynumlayers = numlayers;
    
    myextrudedirection = extrudedirection;

    mybaseshape = innerrawshape;

    mesh();
}

std::shared_ptr<rawshape> rawextrusion::duplicate(void)
{ 
    std::shared_ptr<rawextrusion> out(new rawextrusion);
    *out = *this;

    out->sons = geotools::duplicate(sons);
    out->mybaseshape = mybaseshape->duplicate();

    out->replicatelinks(shared_from_this());

    return out;    
}

void rawextrusion::setphysicalregion(int physreg)
{
    myphysicalregion = physreg;
}

int rawextrusion::getdimension(void)
{
    return mybaseshape->getdimension() + 1;
}

std::string rawextrusion::getname(void)
{
    return "extrusion";
}


std::vector<std::shared_ptr<rawshape>> rawextrusion::getsons(void) 
{ 
    return sons; 
}

std::vector<std::shared_ptr<rawshape>> rawextrusion::getsubshapes(void)
{
    return sons;
}

void rawextrusion::setsubshapes(std::vector<std::shared_ptr<rawshape>> subshapes)
{
    sons = subshapes;
}

int rawextrusion::getphysicalregion(void) 
{ 
    return myphysicalregion; 
}

std::vector<double>* rawextrusion::getcoords(void) 
{ 
    return &mycoords; 
}

std::vector<std::vector<int>>* rawextrusion::getelems(void) 
{ 
    return &myelems; 
}

std::shared_ptr<rawshape> rawextrusion::getpointer(void) 
{ 
    return shared_from_this(); 
}


void rawextrusion::mesh(void)
{
    // The extrusion of a 0D shape (point) is a line:
    if (mybaseshape->getdimension() == 0)
    {
        // Create the point at the top side of the extrusion:
        std::shared_ptr<rawshape> toppoint = mybaseshape->duplicate();
        toppoint->setphysicalregion(-1);
        toppoint->shift(myheight*myextrudedirection[0], myheight*myextrudedirection[1], myheight*myextrudedirection[2]);

        sons = {mybaseshape, toppoint};

        
        int numnodesperlayer = mybaseshape->getcoords()->size()/3;
        
        mycoords.resize(3*numnodesperlayer*mynumlayers);
        // Extruded points are lines (element number 1):
        myelems[1].resize(2*numnodesperlayer*(mynumlayers-1));


        // Place the nodes.
        int index = 0;
        
        std::vector<double>* currentcoords = mybaseshape->getcoords();

        for (int l = 0; l < mynumlayers; l++)
        {
            for (int j = 0; j < numnodesperlayer; j++)
            {
                mycoords[3*index+0] = currentcoords->at(3*j+0) + l/(mynumlayers-1.0)*myheight*myextrudedirection[0];
                mycoords[3*index+1] = currentcoords->at(3*j+1) + l/(mynumlayers-1.0)*myheight*myextrudedirection[1];
                mycoords[3*index+2] = currentcoords->at(3*j+2) + l/(mynumlayers-1.0)*myheight*myextrudedirection[2];
                index++;
            }
        }


        // Place the elements:
        int lineindex = 0;

        std::vector<std::vector<int>>* currentelems = mybaseshape->getelems();

        for (int l = 0; l < mynumlayers-1; l++)
        {
            // Extrude all points:
            for (int p = 0; p < numnodesperlayer; p++)
            {
                myelems[1][2*lineindex+0] = currentelems->at(0)[p] + l*numnodesperlayer;
                myelems[1][2*lineindex+1] = currentelems->at(0)[p] + (l+1)*numnodesperlayer;
                lineindex++;
            }
        }
    }

    // The extrusion of a 1D shape is a quadrangle:
    if (mybaseshape->getdimension() == 1)
    {
        // Create the line at the top side of the extrusion:
        std::shared_ptr<rawshape> topline = mybaseshape->duplicate();
        topline->setphysicalregion(-1);
        topline->shift(myheight*myextrudedirection[0], myheight*myextrudedirection[1], myheight*myextrudedirection[2]);

        // Create the lines at both sides of the quadrangle:
        if ((mybaseshape->getsons()).size() == 2)
        {
            std::shared_ptr<rawshape> leftline = std::shared_ptr<rawline>(new rawline(-1, {mybaseshape->getsons()[0], topline->getsons()[0]}, mynumlayers));
            std::shared_ptr<rawshape> rightline = std::shared_ptr<rawline>(new rawline(-1, {mybaseshape->getsons()[1], topline->getsons()[1]}, mynumlayers));

            sons = {mybaseshape, rightline, topline, leftline};
        }
        else
            sons = {mybaseshape, topline};
            

        int numnodesperlayer = mybaseshape->getcoords()->size()/3;
        int numlines = mybaseshape->getelems()->at(1).size()/2;
        
        mycoords.resize(3*numnodesperlayer*mynumlayers);
        // Extruded lines are quadrangles (element number 3):
        myelems[3].resize(4*numlines*(mynumlayers-1));


        // Place the nodes.
        int index = 0;
        
        std::vector<double>* currentcoords = mybaseshape->getcoords();

        for (int l = 0; l < mynumlayers; l++)
        {
            for (int j = 0; j < numnodesperlayer; j++)
            {
                mycoords[3*index+0] = currentcoords->at(3*j+0) + l/(mynumlayers-1.0)*myheight*myextrudedirection[0];
                mycoords[3*index+1] = currentcoords->at(3*j+1) + l/(mynumlayers-1.0)*myheight*myextrudedirection[1];
                mycoords[3*index+2] = currentcoords->at(3*j+2) + l/(mynumlayers-1.0)*myheight*myextrudedirection[2];
                index++;
            }
        }


        // Place the elements:
        int quadindex = 0;

        std::vector<std::vector<int>>* currentelems = mybaseshape->getelems();

        for (int l = 0; l < mynumlayers-1; l++)
        {
            // Extrude all lines:
            for (int p = 0; p < numlines; p++)
            {
                myelems[3][4*quadindex+0] = currentelems->at(1)[2*p+0] + l*numnodesperlayer;
                myelems[3][4*quadindex+1] = currentelems->at(1)[2*p+1] + l*numnodesperlayer;
                myelems[3][4*quadindex+2] = currentelems->at(1)[2*p+1] + (l+1)*numnodesperlayer;
                myelems[3][4*quadindex+3] = currentelems->at(1)[2*p+0] + (l+1)*numnodesperlayer;
                quadindex++;
            }
        }
    }

    // The extrusion of a 2D shape is a volume:
    if (mybaseshape->getdimension() == 2)
    {
        // Create the face at the top side of the extrusion:
        std::shared_ptr<rawshape> topface = mybaseshape->duplicate();
        topface->setphysicalregion(-1);
        topface->shift(myheight*myextrudedirection[0], myheight*myextrudedirection[1], myheight*myextrudedirection[2]);

        // Get the contour regions:
        std::vector<std::shared_ptr<rawshape>> contourlines = mybaseshape->getsons();
        std::vector<std::shared_ptr<rawshape>> topcontourlines = topface->getsons();

        contourlines = geotools::orient(contourlines);
        topcontourlines = geotools::orient(topcontourlines);

        // Create the contour faces:
        std::vector<std::shared_ptr<rawshape>> contourfaces(contourlines.size());

        std::vector<std::shared_ptr<rawshape>> verticallines(contourlines.size());
        for (int i = 0; i < contourlines.size(); i++)
            verticallines[i] = std::shared_ptr<rawline>(new rawline(-1, {contourlines[i]->getsons()[0], topcontourlines[i]->getsons()[0]}, mynumlayers));
        
        for (int i = 0; i < contourfaces.size(); i++)
            contourfaces[i] = std::shared_ptr<rawquadrangle>(new rawquadrangle(-1, {contourlines[i], verticallines[(i+1)%verticallines.size()], topcontourlines[i], verticallines[i]}));

        sons = geotools::concatenate({{mybaseshape}, contourfaces, {topface}});


        // Get total number of nodes, triangles and quadrangles in the unextruded mesh:
        int numnodesperlayer = mybaseshape->getcoords()->size()/3;
        int numtriangles = mybaseshape->getelems()->at(2).size()/3;
        int numquadrangles = mybaseshape->getelems()->at(3).size()/4;

        mycoords.resize(3*numnodesperlayer*mynumlayers);
        // Extruded triangles are prisms (element number 6):
        myelems[6].resize(6*numtriangles*(mynumlayers-1));
        // Extruded quadrangles are hexahedra (element number 5):
        myelems[5].resize(8*numquadrangles*(mynumlayers-1));


        // Place the nodes.
        int index = 0;
        
        std::vector<double>* currentcoords = mybaseshape->getcoords();

        for (int l = 0; l < mynumlayers; l++)
        {
            for (int j = 0; j < numnodesperlayer; j++)
            {
                mycoords[3*index+0] = currentcoords->at(3*j+0) + l/(mynumlayers-1.0)*myheight*myextrudedirection[0];
                mycoords[3*index+1] = currentcoords->at(3*j+1) + l/(mynumlayers-1.0)*myheight*myextrudedirection[1];
                mycoords[3*index+2] = currentcoords->at(3*j+2) + l/(mynumlayers-1.0)*myheight*myextrudedirection[2];
                index++;
            }
        }


        // Place the elements:
        int prismindex = 0, hexindex = 0;

        std::vector<std::vector<int>>* currentelems = mybaseshape->getelems();

        for (int l = 0; l < mynumlayers-1; l++)
        {
            // Extrude all triangles:
            for (int tri = 0; tri < numtriangles; tri++)
            {
                myelems[6][6*prismindex+0] = currentelems->at(2)[3*tri+0] + l*numnodesperlayer;
                myelems[6][6*prismindex+1] = currentelems->at(2)[3*tri+1] + l*numnodesperlayer;
                myelems[6][6*prismindex+2] = currentelems->at(2)[3*tri+2] + l*numnodesperlayer;
                myelems[6][6*prismindex+3] = currentelems->at(2)[3*tri+0] + (l+1)*numnodesperlayer;
                myelems[6][6*prismindex+4] = currentelems->at(2)[3*tri+1] + (l+1)*numnodesperlayer;
                myelems[6][6*prismindex+5] = currentelems->at(2)[3*tri+2] + (l+1)*numnodesperlayer;
                prismindex++;
            }

            // Extrude all quadrangles:
            for (int quad = 0; quad < numquadrangles; quad++)
            {
                myelems[5][8*hexindex+0] = currentelems->at(3)[4*quad+0] + l*numnodesperlayer;
                myelems[5][8*hexindex+1] = currentelems->at(3)[4*quad+1] + l*numnodesperlayer;
                myelems[5][8*hexindex+2] = currentelems->at(3)[4*quad+2] + l*numnodesperlayer;
                myelems[5][8*hexindex+3] = currentelems->at(3)[4*quad+3] + l*numnodesperlayer;
                myelems[5][8*hexindex+4] = currentelems->at(3)[4*quad+0] + (l+1)*numnodesperlayer;
                myelems[5][8*hexindex+5] = currentelems->at(3)[4*quad+1] + (l+1)*numnodesperlayer;
                myelems[5][8*hexindex+6] = currentelems->at(3)[4*quad+2] + (l+1)*numnodesperlayer;
                myelems[5][8*hexindex+7] = currentelems->at(3)[4*quad+3] + (l+1)*numnodesperlayer;
                hexindex++;
            }
        }
    }
}



