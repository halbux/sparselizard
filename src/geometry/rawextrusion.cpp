#include "rawextrusion.h"


rawextrusion::rawextrusion(int physreg, std::vector<std::shared_ptr<rawshape>> contour, std::vector<std::shared_ptr<rawshape>> innerregions, double height, int numlayers)
{
	myphysicalregion = physreg;

	// Sons will be created while meshing

	// The extruded shape is of higher dimension:
	mydimension = innerregions[0]->getdimension() + 1;

	myheight = height;

	mynumlayers = numlayers;

	myunextrudedregions = innerregions;

	mycontourregions = contour;

	mesh();
}

std::shared_ptr<rawshape> rawextrusion::duplicate(void)
{ 
	std::shared_ptr<rawextrusion> out(new rawextrusion);
	*out = *this;

	out->sons = geotools::duplicate(sons);
	out->myunextrudedregions = geotools::duplicate(myunextrudedregions);
	out->mycontourregions = geotools::duplicate(mycontourregions);

	return out;	
}

void rawextrusion::setphysicalregion(int physreg)
{
	myphysicalregion = physreg;
}

int rawextrusion::getdimension(void)
{
	return mydimension;
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
	return geotools::concatenate({sons,myunextrudedregions,mycontourregions});
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
	if (mydimension-1 == 0)
	{
		// Create the point at the top side of the extrusion:
		std::shared_ptr<rawshape> toppoint = myunextrudedregions[0]->duplicate();
		toppoint->setphysicalregion(-1);
		toppoint->shift(0, 0, myheight);

		sons = {myunextrudedregions[0], toppoint};

		// Create this temporary object to get the mesh:
		std::shared_ptr<rawshape> rawshapeptr = std::shared_ptr<rawline>(new rawline(myphysicalregion, {myunextrudedregions[0], toppoint}, mynumlayers));

		mycoords = *(rawshapeptr->getcoords());
		myelems = *(rawshapeptr->getelems());
	}

	// The extrusion of a 1D shape is a quadrangle:
	if (mydimension-1 == 1)
	{
		myunextrudedregions = geotools::orient(myunextrudedregions);

		// Create the line at the top side of the extrusion:
		std::vector<std::shared_ptr<rawshape>> toplines = geotools::duplicate(myunextrudedregions);
		for (int i = 0; i < toplines.size(); i++)	
		{	
			toplines[i]->setphysicalregion(-1);
			toplines[i]->shift(0, 0, myheight);
		}

		// Create the lines at both sides of the quadrangle (if any):
		if (mycontourregions.size() > 0)
		{
			std::vector<std::shared_ptr<rawshape>> baselinepts = {mycontourregions[0], mycontourregions[1]};
			std::vector<std::shared_ptr<rawshape>> toplinespts = {toplines[0]->getsons()[0], toplines[toplines.size()-1]->getsons()[1]};

			std::shared_ptr<rawshape> leftline = std::shared_ptr<rawline>(new rawline(-1, {baselinepts[0], toplinespts[0]}, mynumlayers));
			std::shared_ptr<rawshape> rightline = std::shared_ptr<rawline>(new rawline(-1, {baselinepts[1], toplinespts[1]}, mynumlayers));

			sons = geotools::orient(geotools::concatenate({myunextrudedregions, {rightline}, geotools::flip(toplines), {leftline}}));
		}
		else
			sons = geotools::orient(geotools::concatenate({myunextrudedregions, geotools::flip(toplines)}));


		///// Now mesh the extruded quadrangle:		

		// Count total number of nodes and lines in the unextruded mesh:
		int numnodes = 0, numlines = 0;

		for (int i = 0; i < myunextrudedregions.size(); i++)
		{
			numnodes += myunextrudedregions[i]->getcoords()->size()/3;
			numlines += myunextrudedregions[i]->getelems()->at(1).size()/2;
		}

		mycoords.resize(3*numnodes*mynumlayers);
		// Extruded lines are quadrangles (element number 3):
		myelems[3].resize(4*numlines*(mynumlayers-1));


		// Place the nodes.
		int index = 0;
		for (int l = 0; l < mynumlayers; l++)
		{
			for (int i = 0; i < myunextrudedregions.size(); i++)
			{	
				std::vector<double>* currentcoords = myunextrudedregions[i]->getcoords();

				for (int j = 0; j < currentcoords->size()/3; j++)
				{
					mycoords[3*index+0] = currentcoords->at(3*j+0);
					mycoords[3*index+1] = currentcoords->at(3*j+1);
					mycoords[3*index+2] = currentcoords->at(3*j+2) + l/(mynumlayers-1.0)*myheight;
					index++;
				}
			}
		}

		int numnodesperlayer = mycoords.size()/3/mynumlayers;

		// Place the elements:
		int lineindex = 0;
		// The node numbers need to be shifted:
		int nodeshift = 0;
		for (int i = 0; i < myunextrudedregions.size(); i++)
		{	
			std::vector<std::vector<int>>* currentelems = myunextrudedregions[i]->getelems();

			// Extrude all lines:
			for (int lin = 0; lin < currentelems->at(1).size()/2; lin++)
			{
				for (int l = 0; l < mynumlayers-1; l++)
				{
					myelems[3][4*lineindex+0] = currentelems->at(1)[2*lin+0] + nodeshift + l*numnodesperlayer;
					myelems[3][4*lineindex+1] = currentelems->at(1)[2*lin+1] + nodeshift + l*numnodesperlayer;
					myelems[3][4*lineindex+2] = currentelems->at(1)[2*lin+1] + nodeshift + (l+1)*numnodesperlayer;
					myelems[3][4*lineindex+3] = currentelems->at(1)[2*lin+0] + nodeshift + (l+1)*numnodesperlayer;
					lineindex++;
				}
			}
			nodeshift += myunextrudedregions[i]->getcoords()->size()/3;
		}

	}

	// The extrusion of a 2D shape is a volume:
	if (mydimension-1 == 2)
	{
		// Create the face at the top side of the extrusion:
		std::vector<std::shared_ptr<rawshape>> topfaces = geotools::duplicate(myunextrudedregions);
		for (int i = 0; i < topfaces.size(); i++)	
		{	
			topfaces[i]->setphysicalregion(-1);
			topfaces[i]->shift(0, 0, myheight);
		}

		// Create the contour faces:
		std::vector<std::shared_ptr<rawshape>> contourfaces(mycontourregions.size());
		for (int i = 0; i < mycontourregions.size(); i++)
			contourfaces[i] = std::shared_ptr<rawshape>(new rawextrusion(-1, {mycontourregions[i]->getsons()[0],mycontourregions[i]->getsons()[1]}, {mycontourregions[i]}, myheight, mynumlayers));

		sons = geotools::concatenate({myunextrudedregions, contourfaces, topfaces});


		// Count total number of nodes, triangles and quadrangles in the unextruded mesh:
		int numnodes = 0, numtriangles = 0, numquadrangles = 0;

		for (int i = 0; i < myunextrudedregions.size(); i++)
		{
			numnodes += myunextrudedregions[i]->getcoords()->size()/3;
			numtriangles += myunextrudedregions[i]->getelems()->at(2).size()/3;
			numquadrangles += myunextrudedregions[i]->getelems()->at(3).size()/4;
		}

		mycoords.resize(3*numnodes*mynumlayers);
		// Extruded triangles are prisms (element number 6):
		myelems[6].resize(6*numtriangles*(mynumlayers-1));
		// Extruded quadrangles are hexahedra (element number 5):
		myelems[5].resize(8*numquadrangles*(mynumlayers-1));


		// Place the nodes.
		int index = 0;
		for (int l = 0; l < mynumlayers; l++)
		{
			for (int i = 0; i < myunextrudedregions.size(); i++)
			{	
				std::vector<double>* currentcoords = myunextrudedregions[i]->getcoords();

				for (int j = 0; j < currentcoords->size()/3; j++)
				{
					mycoords[3*index+0] = currentcoords->at(3*j+0);
					mycoords[3*index+1] = currentcoords->at(3*j+1);
					mycoords[3*index+2] = currentcoords->at(3*j+2) + l/(mynumlayers-1.0)*myheight;
					index++;
				}
			}
		}

		int numnodesperlayer = mycoords.size()/3/mynumlayers;

		// Place the elements:
		int prismindex = 0, hexindex = 0;
		// The node numbers need to be shifted:
		int nodeshift = 0;
		for (int i = 0; i < myunextrudedregions.size(); i++)
		{	
			std::vector<std::vector<int>>* currentelems = myunextrudedregions[i]->getelems();

			// Extrude all triangles:
			for (int tri = 0; tri < currentelems->at(2).size()/3; tri++)
			{
				for (int l = 0; l < mynumlayers-1; l++)
				{
					myelems[6][6*prismindex+0] = currentelems->at(2)[3*tri+0] + nodeshift + l*numnodesperlayer;
					myelems[6][6*prismindex+1] = currentelems->at(2)[3*tri+1] + nodeshift + l*numnodesperlayer;
					myelems[6][6*prismindex+2] = currentelems->at(2)[3*tri+2] + nodeshift + l*numnodesperlayer;
					myelems[6][6*prismindex+3] = currentelems->at(2)[3*tri+0] + nodeshift + (l+1)*numnodesperlayer;
					myelems[6][6*prismindex+4] = currentelems->at(2)[3*tri+1] + nodeshift + (l+1)*numnodesperlayer;
					myelems[6][6*prismindex+5] = currentelems->at(2)[3*tri+2] + nodeshift + (l+1)*numnodesperlayer;
					prismindex++;
				}
			}

			// Extrude all quadrangles:
			for (int quad = 0; quad < currentelems->at(3).size()/4; quad++)
			{
				for (int l = 0; l < mynumlayers-1; l++)
				{
					myelems[5][8*hexindex+0] = currentelems->at(3)[4*quad+0] + nodeshift + l*numnodesperlayer;
					myelems[5][8*hexindex+1] = currentelems->at(3)[4*quad+1] + nodeshift + l*numnodesperlayer;
					myelems[5][8*hexindex+2] = currentelems->at(3)[4*quad+2] + nodeshift + l*numnodesperlayer;
					myelems[5][8*hexindex+3] = currentelems->at(3)[4*quad+3] + nodeshift + l*numnodesperlayer;
					myelems[5][8*hexindex+4] = currentelems->at(3)[4*quad+0] + nodeshift + (l+1)*numnodesperlayer;
					myelems[5][8*hexindex+5] = currentelems->at(3)[4*quad+1] + nodeshift + (l+1)*numnodesperlayer;
					myelems[5][8*hexindex+6] = currentelems->at(3)[4*quad+2] + nodeshift + (l+1)*numnodesperlayer;
					myelems[5][8*hexindex+7] = currentelems->at(3)[4*quad+3] + nodeshift + (l+1)*numnodesperlayer;
					hexindex++;
				}
			}
			nodeshift += myunextrudedregions[i]->getcoords()->size()/3;
		}
	}
}



