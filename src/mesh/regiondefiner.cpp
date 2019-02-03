#include "regiondefiner.h"


void regiondefiner::defineskinregions(void)
{
	myphysicalregions->errorundefined(toskin);

	std::vector<double>* nodecoords = mynodes->getcoordinates();

	// Loop on all skin requests:
	for (int i = 0; i < toskin.size(); i++)
	{
		physicalregion* newphysreg = myphysicalregions->get(skins[i]);
		physicalregion* curphysreg = myphysicalregions->get(toskin[i]);
		
		if (curphysreg->getelementdimension() == 0)
		{
			std::cout << "Error in 'regiondefiner' object: cannot get the skin of point elements" << std::endl;
			abort();
		}
		
		int physregdim = curphysreg->getelementdimension();
		// The skin has a one lower dimension:
		int skindim = physregdim-1;
		
		std::vector<std::vector<int>>* curelems = curphysreg->getelementlist();
		
		// Loop on all skin element types:
		for (int skinelemtype = 0; skinelemtype <= 7; skinelemtype++)
		{
			element myskinelement(skinelemtype);
			if (myskinelement.getelementdimension() != skindim)
				continue;

			// Vector to count the number of times a skin element appears in an element:
			std::vector<int> numoccurences(myelements->count(skinelemtype),0);
				
			// Loop on all element types in the physical region:
			for (int elemtype = 0; elemtype <= 7; elemtype++)
			{
				std::vector<int>* curelemtype = &(curelems->at(elemtype));
			
				int numelems = curelemtype->size();
				if (numelems == 0)
					continue;
		
				element myelement(elemtype);
				
				int numsubtypeelems = myelement.counttype(skinelemtype);
				
				// Loop on all elements:
				for (int elem = 0; elem < numelems; elem++)
				{
					int curelem = curelemtype->at(elem);
				
					for (int subelem = 0; subelem < numsubtypeelems; subelem++)
						numoccurences[myelements->getsubelement(skinelemtype, elemtype, curelem, subelem)]++;
				}
			}
			
			// All single occurences are skin elements:
			for (int j = 0; j < numoccurences.size(); j++)
			{
				if (numoccurences[j] == 1)
					newphysreg->addelement(skinelemtype, j);
			}
		}
		newphysreg->removeduplicatedelements();
	}
}

void regiondefiner::defineboxregions(void)
{
	myphysicalregions->errorundefined(tobox);
	
	std::vector<double>* nodecoords = mynodes->getcoordinates();

	// Loop on all box selection requests:
	for (int i = 0; i < tobox.size(); i++)
	{
		std::vector<double> boxlimit = boxlimits[i];
		
		physicalregion* newphysreg = myphysicalregions->get(boxed[i]);
		physicalregion* curphysreg = myphysicalregions->get(tobox[i]);
		
		std::vector<std::vector<int>>* curelems = curphysreg->getelementlist();
		
		// Loop on all element types:
		for (int elemtype = 0; elemtype <= 7; elemtype++)
		{
			if (curelems->at(elemtype).size() == 0)
				continue;
				
			element myelement(elemtype);
		
			std::vector<int>* curelemtype = &(curelems->at(elemtype));
		
			// Loop on all subelement types:
			for (int subelemtype = 0; subelemtype <= 7; subelemtype++)
			{
				// Skip all subelements that are not of the requested dimension:
				element mysubelement(subelemtype);
				int numnodes = mysubelement.countcurvednodes();
				int numsubtypeelems = myelement.counttype(subelemtype);

				if (mysubelement.getelementdimension() != boxelemdims[i] || numsubtypeelems == 0)
					continue;
			
				// Loop on all elements:
				for (int elem = 0; elem < curelemtype->size(); elem++)
				{
					int curelem = curelemtype->at(elem);

					for (int subelem = 0; subelem < numsubtypeelems; subelem++)
					{
						int cursubelem = myelements->getsubelement(subelemtype, elemtype, curelem, subelem);
				
						// Check if the coordinates of all nodes in the subelement are within the box limits:
						bool isinlimits = true;
						for (int node = 0; node < numnodes; node++)
						{
							int curnode = myelements->getsubelement(0, subelemtype, cursubelem, node);
				
							double curnodex = nodecoords->at(3*curnode+0);
							double curnodey = nodecoords->at(3*curnode+1);
							double curnodez = nodecoords->at(3*curnode+2);
				
							if (curnodex < boxlimit[0] || curnodex > boxlimit[1] || curnodey < boxlimit[2] || curnodey > boxlimit[3] || curnodez < boxlimit[4] || curnodez > boxlimit[5])
							{
								isinlimits = false;
								break;
							}
						}
						if (isinlimits)
							newphysreg->addelement(subelemtype, cursubelem);
					}
				}
			}
		}
		newphysreg->removeduplicatedelements();
	}
}

regiondefiner::regiondefiner(nodes& inputnodes, elements& inputelems, physicalregions& inputphysregs)
{
	mynodes = &inputnodes;
	myelements = &inputelems;
	myphysicalregions = &inputphysregs;
}

void regiondefiner::regionskin(int newphysreg, int physregtoskin)
{	
	skins.push_back(newphysreg);
	toskin.push_back(physregtoskin);
}

void regiondefiner::boxselection(int newphysreg, int selecteddim, std::vector<double> boxlimit, int physregtobox)
{
	if (boxlimit.size() != 6)
	{
		std::cout << "Error in 'regiondefiner' object: expected a vector of length 6 for the box limits {x1,x2,y1,y2,z1,z2}" << std::endl;
		abort();
	}
	if (selecteddim > 3 || selecteddim < 0)
	{
		std::cout << "Error in 'regiondefiner' object: dimension of the elements to select cannot be " << selecteddim << std::endl;
		abort();
	}
	
	boxed.push_back(newphysreg);
	tobox.push_back(physregtobox);
	boxelemdims.push_back(selecteddim);
	
	// Make the box limit slightly larger to remove the roundoff noise issues:
	boxlimit[0] -= boxlimit[0]*roundoffnoise; boxlimit[2] -= boxlimit[2]*roundoffnoise; boxlimit[4] -= boxlimit[4]*roundoffnoise;
 	boxlimit[1] += boxlimit[1]*roundoffnoise; boxlimit[3] += boxlimit[3]*roundoffnoise; boxlimit[5] += boxlimit[5]*roundoffnoise;
 	
	boxlimits.push_back(boxlimit);
}


void regiondefiner::defineregions(void)
{
	defineskinregions();
	defineboxregions();
}


