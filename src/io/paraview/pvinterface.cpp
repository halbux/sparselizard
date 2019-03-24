#include "pvinterface.h"


void pvinterface::writetofile(std::string name, iodata datatowrite)
{
    // Get the file name without the .vtk extension:
    std::string namenoext = name.substr(0, name.size()-4);
    
    mystring myname(name);
    std::string viewname = myname.getstringwhileletter();
    
	// 'file' cannot take a std::string argument --> name.c_str():
	std::ofstream outfile (name.c_str());
	if (outfile.is_open())
	{
		// To write all doubles with enough digits to the file:
		outfile << std::setprecision(17);
		
		// Write the header:
		outfile << "# vtk DataFile Version 4.2\n";
		outfile << namenoext+"\n";
		outfile << "ASCII\n";
		outfile << "DATASET UNSTRUCTURED_GRID\n\n";
		
		// Write the points section.
		int numnodes = datatowrite.countcoordnodes();
		outfile << "POINTS " << numnodes << " double\n";
		for (int i = 0; i < 8; i++)
		{
			if (datatowrite.ispopulated(i) == false)
				continue;
				
			std::vector<densematrix> curcoords = datatowrite.getcoordinates(i);
			double* xvals = curcoords[0].getvalues();
			double* yvals = curcoords[1].getvalues();
			double* zvals = curcoords[2].getvalues();
			
			int index = 0;
			for (int elem = 0; elem < curcoords[0].countrows(); elem++)
			{
				for (int node = 0; node < curcoords[0].countcolumns(); node++)
				{
					outfile << xvals[index] << " " << yvals[index] << " " << zvals[index] << " ";
					index++;
				}
				outfile << "\n";
			}
		}
		outfile << "\n";
		
		// Write the cells section:
		int numelems = datatowrite.countelements();
		outfile << "CELLS " << numelems << " " << numnodes + numelems << "\n";
		
		int nodenum = 0;
		for (int i = 0; i < 8; i++)
		{
			if (datatowrite.ispopulated(i) == false)
				continue;
				
			std::vector<densematrix> curcoords = datatowrite.getcoordinates(i);
			
			for (int elem = 0; elem < curcoords[0].countrows(); elem++)
			{
				outfile << curcoords[0].countcolumns() << " ";
				for (int node = 0; node < curcoords[0].countcolumns(); node++)
				{
					outfile << nodenum << " ";
					nodenum++;
				}
				outfile << "\n";
			}
		}
		outfile << "\n";
		
		// Write the cell types section:
		outfile << "CELL_TYPES " << numelems << "\n";
		for (int i = 0; i < 8; i++)
		{
			if (datatowrite.ispopulated(i) == false)
				continue;
				
			element myelem(i, datatowrite.getinterpolorder());
				
			std::vector<densematrix> curcoords = datatowrite.getcoordinates(i);
			
			for (int elem = 0; elem < curcoords[0].countrows(); elem++)
				outfile << converttoparaviewelementtypenumber(myelem.getcurvedtypenumber()) << "\n";
		}
		outfile << "\n";
		
		// Write the data section:
		outfile << "POINT_DATA " << numnodes << "\n";
		outfile << "\n";
		
		// Write the scalar data section (if any):
		if (datatowrite.isscalar() == true)
		{
			outfile << "SCALARS " << viewname << " double\n";
			outfile << "LOOKUP_TABLE default" << "\n";
			for (int i = 0; i < 8; i++)
			{
				if (datatowrite.ispopulated(i) == false)
					continue;
					
				densematrix scaldat = datatowrite.getdata(i)[0];
				double* scalvals = scaldat.getvalues();
				
				for (int i = 0; i < scaldat.count(); i++)
					outfile << scalvals[i] << "\n";
			}
			outfile << "\n";
		}
		
		// Write the vector data section (if any):
		if (datatowrite.isscalar() == false)
		{
			outfile << "VECTORS " << viewname << " double\n";
			for (int i = 0; i < 8; i++)
			{
				if (datatowrite.ispopulated(i) == false)
					continue;
					
				std::vector<densematrix> vecdat = datatowrite.getdata(i);
				double* compxvals = vecdat[0].getvalues();
				double* compyvals = vecdat[1].getvalues();
				double* compzvals = vecdat[2].getvalues();
				
				for (int i = 0; i < vecdat[0].count(); i++)
					outfile << compxvals[i] << " " << compyvals[i] << " " << compzvals[i] << "\n";
			}
			outfile << "\n";
		}
		
		outfile.close();
	}
	else 
	{
		std::cout << "Unable to write to file " << name << " or file not found" << std::endl;
		abort();
	}
}

int pvinterface::converttoparaviewelementtypenumber(int ourtypenumber)
{
	switch (ourtypenumber)
	{
		// Point:
		case 0:
			return 1;
			
		// Line order 1:
		case 1:
			return 3;
		// Triangle order 1:
		case 2:
			return 5;
		// Quadrangle order 1:
		case 3:
			return 9;
		// Tetrahedron order 1:
		case 4:
			return 10;
		// Hexahedron order 1:
		case 5:
			return 12;
		// Prism order 1:
		case 6:
			return 13;
		// Pyramid order 1:
		case 7:
			return 14;
			
		// Line order 2:
		case 8:
			return 21;
		// Triangle order 2:
		case 9:
			return 22;
		// Quadrangle order 2:
		case 10:
			return 23;
		// Tetrahedron order 2:
		case 11:
			return 24;
		// Hexahedron order 2:
		case 12:
			return 25;
		// Prism order 2:
		case 13:
			return 26;
		// Pyramid order 2:
		case 14:
			return 27;
            
		default:
			std::cout << "Error in 'pvinterface' namespace: trying to use a ParaView element that is undefined in this code." << std::endl;
			abort();
	}
}

