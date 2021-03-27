// This code shows how to perfom a magnetostatic analysis on a 2D cross-section of a
// rotating PMSM (permanent magnet synchronous) electric motor. The motor torque is 
// calculated based on the virtual work principle for an increasing mechanical angle.
//
// Antiperiodicity is used to reduce the computational domain to only 45 degrees of 
// the total geometry (the motor has 4 pairs of poles). In order to link the rotor and
// stator domains at their interface a general mortar-based continuity condition is used.
// This allows to work with the non-matching mesh at the interface when the rotor moves.  

#include <fstream>
#include <experimental/filesystem>
#include <iostream>
#include "sparselizard.h"


using namespace std;
using namespace sl;

expression normalize(expression);
expression polar(expression);
expression polarscalar(expression);

expression projection(expression e, expression a);

void createPath(std::string path) {
    std::string command = "mkdir -p " + path;
    bool success = system(command.c_str());
    // namespace fs = std::experimental::filesystem; // In C++17 use std::filesystem.

    // std::error_code ec;
    // bool success = fs::create_directories(path, ec);

    if (!success) {
        std::cout << "Couldn't create directory path: " << path << std::endl; // Fun fact: In case of success ec.message() returns "The operation completed successfully." using vc++.
    }
}

void writeFile(std::string filename, std::string data) {
    ofstream myfile;
    myfile.open (filename.c_str());
    myfile << data;
    myfile.close();
}

struct DataPack {
    expression A;
    expression B;
    expression F;
    double T;
    double TBottom;
    double TTop;
    double TIron;
    double TBM;
    double TTM;
    double alpha;
    double ironAlpha;

    std::string name;
    int allRegion;

    DataPack(std::string name) {
        this->name = name;
    }

    void saveData() {
        cout << "Saving datapack to folder " << this->name << endl;

        std::string path = "results/" + this->name;
        createPath(path);

        createPath(path + "/B");
        createPath(path + "/A");
        createPath(path + "/F");
        createPath(path + "/T");

        int all = this->allRegion;
        int alpha = this->alpha;
        int ironAlpha = this->ironAlpha;
        
        this->A.write(all, path + "/A/AField"+std::to_string((int)alpha)+".vtu", 2);
        this->B.write(all, path + "/B/BField"+std::to_string((int)alpha)+std::to_string((int)ironAlpha)+".vtu", 2);
        this->F.write(all, path + "/F/Force"+std::to_string((int)alpha)+".vtu", 2);
        std::string data = "Bottom,Iron,Top,Bottom Magnets,Top Magnets\r\n" +
            std::to_string(this->TBottom) + "," + std::to_string(this->TIron) + "," + std::to_string(this->TTop) + "," + std::to_string(this->TBM) + "," + std::to_string(this->TTM);
        writeFile(path + "/T/torgues"+std::to_string((int)alpha)+std::to_string((int)ironAlpha)+".csv", data);
        // polar(this->B).write(all, path + "TestPolarB"+std::to_string((int)alpha)+".vtu", 2);
        // polar(this->F).write(all, path + "TestPolarF"+std::to_string((int)alpha)+".vtu", 2);
        // polarscalar(this->B).write(all, path + "TestPolarScalarB"+std::to_string((int)alpha)+".vtu", 2);
    }
};

void writeVec(std::vector<double> v, std::string filename) {
    ofstream myfile(filename);
    int vsize = v.size();
    for (int n=0; n<vsize; n++)
    {
        myfile << v[n] << endl;
    }
}

// Input is rotor angular position in degrees. Output is torque in Nm.
DataPack* sparselizard(mesh mymesh, double alpha = 0, double ironAlpha = 0, bool autoload = false, std::string name="December")
{   
    // Magnet Power
    double MagnetB = 1;

    int stator = 1, firstmagnetup = 2, firstmagnetdown = 3,
     ironbars = 4,
     secondmagnetup = 6, secondmagnetdown = 7,
    //  secondmagnetair = 8, ironair = 9,
     topStator = 20;
    // int IronTop = 21, SecondMagnetBottom = 22;
    // int topStatorTop = 24, bottomStatorBottom = 23;
    
    // Define new physical regions for convenience:
    int magnetsdown = selectunion({secondmagnetdown, firstmagnetdown});
    int magnetsup = selectunion({secondmagnetup, firstmagnetup});
    int magnets = selectunion({firstmagnetdown, firstmagnetup, secondmagnetdown, secondmagnetup});
    int iron = selectunion({ironbars, stator, topStator});
    int topPart = selectunion({ secondmagnetup, secondmagnetdown, topStator });
    int bottomPart = selectunion({ firstmagnetup, firstmagnetdown, stator, ironbars });
    int rotor = selectunion({secondmagnetup, secondmagnetdown});
    int all = selectall();

    int bottommagnets = selectunion({firstmagnetup, firstmagnetdown});
    int topmagnets = selectunion({secondmagnetup, secondmagnetdown});
    int middlesection = selectunion({ironbars});

    // Define a spanning tree to gauge the magnetic vector potential (otherwise the matrix is singular).
    // Start growing the tree from the regions with constrained potential vector (here the contour): 
    spanningtree spantree({topStator});
    // Write it for illustration:
    // spantree.write("results/spantree.pos");

    // Nodal shape functions 'h1' for the z component of the vector potential.
    field az("hcurl", spantree),
         x("x"), y("y"), z("z");


    // maybe  ? 
    az.setgauge(all);
    // Use interpolation order 2:
    az.setorder(all, 1);
    // Put a magnetic wall
    // az.setconstraint(topStator);
    // az.setconstraint(bottomStatorBottom);

    // The remanent induction field in the magnet is 0.5 Tesla perpendicular to the magnet:
    expression normedradialdirection = array3x1(0,0,1);
    expression bremanent = MagnetB * normedradialdirection;
    expression minusnormedradialdirection = array3x1(0,0,-1);
    expression minusbremanent = MagnetB * minusnormedradialdirection;

    // Vacuum magnetic permeability [H/m]: 
    double mu0 = 4.0*getpi()*1e-7;
    // Define the permeability in all regions.
    //
    // Taking into account saturation and measured B-H curves can be easily done
    // by defining an expression based on a 'spline' object (see documentation).
    //
    parameter mu;
    mu|all = mu0;
    // Overwrite on non-magnetic regions:
    mu|magnets = mu0;
    mu|iron = 2000*mu0;

    formulation magnetostatics;
    // The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
    magnetostatics += integral(all, 1/mu* curl(dof(az)) * curl(tf(az)) );

    // Add the remanent magnetization of the rotor magnet:
    magnetostatics += integral(magnetsdown, -1/mu* minusbremanent * curl(tf(az)));
    magnetostatics += integral(magnetsup, -1/mu* bremanent * curl(tf(az)));
    // magnetostatics += continuitycondition(IronTop, SecondMagnetBottom, az, az, 1, false);
    if(autoload == true) {
        az.loadraw("cache/aaz"+std::to_string((int)alpha)+std::to_string((int)ironAlpha)+".slz.gz", true);
    } else {
        std::cout << "Solving..." << std::endl;
        solve(magnetostatics);
        std::cout << "Solved!" << std::endl;
        az.writeraw(all, "cache/aaz"+std::to_string((int)alpha)+std::to_string((int)ironAlpha)+".slz.gz", true);
    }

    // phiProjection.setdata(all, )


    // expression polito = polarscalar(curl(az));
    // double integration = polito.integrate(firstmagnetup, 4);
    // double integrationdown = polito.integrate(firstmagnetdown, 4);

    // polito.write(bottommagnets, "results/bottommagnets"+std::to_string((int)alpha)+".vtu");
    // cout << "First row INTEGRATION  UP : " << polito.integrate(firstmagnetup, 4) << endl;
    // cout << "First row INTEGRATION DOWN: " << polito.integrate(firstmagnetdown, 4) << endl;
    // cout << " ---------------------------------------------------------------- " << integration << endl;
    // cout << " INTEGRATION  IRON : " << polito.integrate(ironbars, 4) << endl;
    // // cout << " INTEGRATION : " << integrationdown << endl;
    // cout << " ---------------------------------------------------------------- " << integration << endl;
    // cout << " Second row INTEGRATION  UP : " << polito.integrate(secondmagnetup, 4) << endl;
    // cout << " Second row INTEGRATION DOWN: " << polito.integrate(secondmagnetdown, 4) << endl;


    // The magnetostatic force acting on the motor is computed below.

    // This field will hold the x and y component of the magnetic forces:
    field magforce("h1xyz");
    magforce.setorder(all, 2);

    // The magnetic force is projected on field 'magforce' on the solid stator region.
    // This is done with a formulation of the type dof*tf - force calculation = 0.
    formulation forceprojection;

    forceprojection += integral(all, dof(magforce)*tf(magforce));
    forceprojection += integral(all, -predefinedmagnetostaticforce(tf(magforce, all), 1/mu*curl(az), mu));

    solve(forceprojection);

    expression magforcescalar = polarscalar((expression)magforce);
    // magforcescalar.write(bottomPart, "results/bottomMagForce"+std::to_string((int)alpha)+".vtu");


    // Calculate the torque:
    expression leverarm = array3x1(x,y,0);
    expression torgue = compz(crossproduct(leverarm, magforce));

    double T = torgue.integrate(bottomPart, 5);
    double TBottom = T;
    double TTop = torgue.integrate(topPart, 5);
    double TIron = torgue.integrate(middlesection, 5);
    double TBM = torgue.integrate(bottommagnets, 5);
    double TTM = torgue.integrate(topmagnets, 5);

    double scaleFactor = 1;
    cout << "TorgueBottom: " << TBottom   * scaleFactor << endl;
    cout << "TorgueIron: " << TIron   * scaleFactor << endl;
    cout << "TorgueTop: " << TTop   * scaleFactor << endl;

    DataPack* pack = new DataPack(name);
    pack->T = T* scaleFactor;
    pack->TBottom = T* scaleFactor;
    pack->TTop = TTop* scaleFactor;
    pack->TIron = TIron* scaleFactor;
    pack->TBM = TBM* scaleFactor;
    pack->TTM = TTM* scaleFactor;

    pack->F = (expression)magforce;
    pack->A = (expression)az;
    pack->B = curl(az);
    pack->alpha = alpha;
    pack->ironAlpha = ironAlpha;
    pack->allRegion = selectall();

    // pack->saveData();

    return pack;
}

std::string getFilename(int alpha,int ironAlpha,int argc,char* argv[]) {
    // FILE MESH 2.2 .msh
    std::string num = ((alpha <= 9) ? "0" : "") + (std::to_string((int)(alpha)));
    std::string numiron = ((ironAlpha <= 9) ? "0" : "") + (std::to_string((int)(ironAlpha)));
    // std::string test = "models/New.msh";
    std::string filename;
    if ( argc < 4 ) {
        // take file from args
        filename = (std::string)"models/" + argv[1];
        std::cout << "Loading model: " << filename << endl;
    } else {
        filename = (std::string)"models/parts/" + argv[1] + num + "A" + numiron + ".msh";
        std::cout << "Loading model(autoname): " << filename << endl;
    }
    return filename;
}

int main(int argc, char *argv[])
{   
    PetscInitialize(0,{},0,0);

    wallclock clk;

    int realIronAngle = 3;
    int realAngle = 5;

    int start = 0;
    int pieces = 0;
    int alphaStep = 5;
    
    int ironstart = 0;
    int ironpieces = 0;
    int ironStep = 5;

    std::string folderName = "December";

    if(argc > 3) {
        alphaStep = atoi(argv[2]);
        ironStep = atoi(argv[3]);
        pieces = atoi(argv[4]);
        ironpieces = atoi(argv[5]);
        folderName = argv[6];
    } else if(argc == 3) {
        folderName = argv[2];
    } else if(argc == 1) {
        std::cout 
        << "Use this software as follows:" << endl
        << "Example command to load only one mesh file" << endl << endl
        << "./loadem.sh FilenameToLoad" << endl << endl
        << "Example command to load only one mesh file save to folderName (default:December)" << endl << endl
        << "./loadem.sh FilenameToLoad FolderName" << endl << endl
        << "Example command to load only multiple mesh file of the name convention: FilenameToLoad{alpha}A{ironAlpha}.msh" << endl << endl << endl
        << "./loadem.sh FilenameToLoad alphaStep ironStep pieces ironpieces" << endl << endl << endl;
        return 0;
    }
    std::cout << "Settings" << endl 
    << "FolderName: " << folderName << endl 
    << "AlphaStep: " << alphaStep << endl 
    << "IronStep: " << ironStep << endl
    << "Pieces: " << pieces << endl
    << "IronPieces: " << ironpieces << endl;

    DataPack* pack;
    // vector<double> torgues;
    // vector<double> torguesBM;
    // vector<double> torguesTM;
    // vector<double> torguesIron;


    // make only one loop 
    for (double ironAlpha = ironstart; ironAlpha <= ironstart+(ironpieces*ironStep); ironAlpha += ironStep){
        for (double alpha = start; alpha <= start+(pieces*alphaStep); alpha += alphaStep)
        {
        // double alpha = 0 - ironAlpha;
        if (alpha < 0) alpha += 18;
        if (alpha < 0) alpha += 18;
            // settime(alpha - start+0.0001);
            std::string filename = getFilename(alpha, ironAlpha, argc, argv);
            mesh mymesh(filename);
            // mymesh.scale(regionall(), 0.01,0.01,0.01);
            pack = sparselizard(mymesh, alpha, ironAlpha, false, folderName);
            pack->saveData();
            
            // torgues.push_back(pack->T);
            // torguesBM.push_back(pack->TBottom);
            // torguesIron.push_back(pack->TIron);
            // torguesTM.push_back(pack->TTop);

            // std::cout << "Iron Angle\t| Mechanical angle [degrees] and torque [Nm]:" << std::endl;
            // std::cout << ironAlpha << "\t\t|" << alpha << "\t\t\t\t\t|" << pack->T << std::endl;   
            
            delete pack;
        }

        // writeVec(torguesBM, "results/IronB"+ std::to_string((int)(ironAlpha))+"A"+ std::to_string((int)(alphaStep))+"torguesTBottom.csv");
        // writeVec(torguesTM, "results/IronB"+ std::to_string((int)(ironAlpha))+"A"+ std::to_string((int)(alphaStep))+"torguesTTop.csv");
        // writeVec(torguesIron, "results/IronB"+ std::to_string((int)(ironAlpha))+"A"+ std::to_string((int)(alphaStep))+"torguesTIron.csv");
        
        // torgues.clear();
        // torguesBM.clear();
        // torguesTM.clear();
        // torguesIron.clear();
    }

    clk.print("Total run time:");

    PetscFinalize();
    
    return 0;
}



expression normalize(expression e) {
    return e / norm(e);
}
expression polar(expression e) {
    field x("x"), y("y"), z("z");
    // expression proj = e * array3x1(-y,x,0);
    expression proj = projection(e, array3x1(-y,x,0));
    return proj;
}

expression polarscalar(expression e) {
    field x("x"), y("y"), z("z");
    expression proj = projection(e, array3x1(-y,x,0));
    return norm(proj);
}

expression projection(expression e, expression a) {
    return ((e * a) / (a * a)) * a;
}
