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
#include <vector>

using namespace std;
using namespace sl;

expression normalize(expression);
expression polar(expression);
expression polarscalar(expression);

expression projection(expression e, expression a);

const char* getVar(int argc, char* argv[], int position) {
    if(argc > position) {
        return argv[position];
    } else {
        cout
            << "Use this software as follows:" << endl
        << "Example command to load only one mesh file" << endl << endl
        << "./loadem.sh FilenameToLoad" << endl << endl
        << "Example command to load only one mesh file save to folderName (default:December)" << endl << endl
        << "./loadem.sh FilenameToLoad FolderName" << endl << endl
        << "Example command to load only multiple mesh file of the name convention: FilenameToLoad{alpha}A{ironAlpha}.msh" << endl << endl << endl
        << "./loadem.sh FilenameToLoad alphaStep ironStep pieces ironpieces" << endl << endl;
        exit(1);
    }
}
void writeVec(vector<double> v, string filename) {
    ofstream myfile(filename);
    int vsize = v.size();
    for (int n=0; n<vsize; n++)
    {
        myfile << v[n] << endl;
    }
}

string twodigits(int num) {
    return ((num <= 9) ? "0" : "") + (to_string((int)(num)));
}



void createPath(string path) {
    string command = "mkdir -p " + path;
    bool success = system(command.c_str());
    // namespace fs = experimental::filesystem; // In C++17 use filesystem.

    // error_code ec;
    // bool success = fs::create_directories(path, ec);

    if (!success) {
        //cout << "Couldn't create directory path: " << path << endl; // Fun fact: In case of success ec.message() returns "The operation completed successfully." using vc++.
    }
}

void writeFile(string filename, string data) {
    ofstream myfile;
    myfile.open (filename.c_str());
    myfile << data;
    myfile.close();
}

struct DataArray {
    vector<double> array;
    string name;

    DataArray(string n) {

        this->name = n;
    }

    void add(double x) {
        this->array.push_back(x);
    }

    void save(string filename) {
        string path = "results/" + this->name;
        createPath(path);
        createPath(path + "/Arrays");
        writeVec(this->array, path + "/Arrays/" + filename + ".csv");
    }
};

struct DataPack {
    expression A;
    expression B;
    expression H;
    expression F;
    double T;
    double TBottom;
    double TTop;
    double TIron;
    double TBM;
    double TTM;
    double alpha;
    double ironAlpha;

    string name;
    int allRegion;

    DataPack(string name) {
        this->name = name;
    }

    void saveData() {
        cout << "Saving datapack to folder " << this->name << endl;

        string path = "results/" + this->name;
        createPath(path);

        createPath(path + "/B");
        createPath(path + "/A");
        createPath(path + "/F");
        createPath(path + "/T");
        createPath(path + "/H");

        int all = this->allRegion;
        int alpha = this->alpha;
        int ironAlpha = this->ironAlpha;
        
        this->A.write(all, path + "/A/AField"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        this->B.write(all, path + "/B/BField"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        this->H.write(all, path + "/H/HField"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        this->F.write(all, path + "/F/Force"+twodigits(ironAlpha)+twodigits(alpha)+".vtu", 2);
        string data = "Bottom,Iron,Top,Bottom Magnets,Top Magnets\r\n" +
            to_string(this->TBottom) + "," + to_string(this->TIron) + "," + to_string(this->TTop) + "," + to_string(this->TBM) + "," + to_string(this->TTM);
        writeFile(path + "/T/torgues"+twodigits(alpha)+twodigits(ironAlpha)+".csv", data);
        // polar(this->B).write(all, path + "TestPolarB"+twodigits(alpha)+".vtu", 2);
        // polar(this->F).write(all, path + "TestPolarF"+twodigits(alpha)+".vtu", 2);
        // polarscalar(this->B).write(all, path + "TestPolarScalarB"+twodigits(alpha)+".vtu", 2);
    }
};


// Input is rotor angular position in degrees. Output is torque in Nm.
DataPack* sparselizard(mesh mymesh, double alpha = 0, double ironAlpha = 0, bool autoload = false, string name="December", bool timeit = false)
{
    wallclock clk;
    // Magnet Power
    double MagnetB = 1;

    int bottomStator = 1,
        firstmagnetup = 2,
        firstmagnetdown = 3,
        ironbars = 4,
        secondmagnetup = 6,
        secondmagnetdown = 7,
        air = 10,
        topStator = 20;

    // Define new physical regions for convenience:
    int magnetsdown = selectunion({secondmagnetdown, firstmagnetdown});
    int magnetsup = selectunion({secondmagnetup, firstmagnetup});
    int magnets = selectunion({firstmagnetdown, firstmagnetup, secondmagnetdown, secondmagnetup});
    int iron = selectunion({ironbars, bottomStator, topStator});
    int topPart = selectunion({ secondmagnetup, secondmagnetdown, topStator });
    int bottomPart = selectunion({ firstmagnetup, firstmagnetdown, bottomStator });
    int rotor = selectunion({secondmagnetup, secondmagnetdown});
    int all = selectall();
    int bottommagnets = selectunion({firstmagnetup, firstmagnetdown});
    int topmagnets = selectunion({secondmagnetup, secondmagnetdown});
    int middlesection = selectunion({ironbars});

cout << "Calculating Spanning Tree" << endl;
    // Define a spanning tree to gauge the magnetic vector potential (otherwise the matrix is singular).
    // Start growing the tree from the regions with constrained potential vector (here the contour): 
//    spanningtree spantree({ bottomStator, topStator });
    // change the span tree to magnet areas
    spanningtree spantree({ magnetsdown, magnetsup });


    // Write it for illustration:
    // spantree.write("results/spantree.pos");
cout << "Setting field" << endl;
    // Nodal shape functions 'h1' for the z component of the vector potential.
    field az("hcurl", spantree),
         x("x"), y("y"), z("z");


    // maybe  ? 
    az.setgauge(all);
    // Use interpolation order 2:
    az.setorder(all, 2);
//    az.setorder(iron, 3);
    // Put a magnetic wall
    // az.setconstraint(topStator);
    // az.setconstraint(bottomStator);

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

// Β-H curve can be performed like below
//std :: vector<double> temperature = {273,300,320,340};
//std :: vector<double> youngsmodulus = {5e9,4e9,2.5e9,1e9};
//spline spl(temperature, youngsmodulus);
//// The spline object can also be created from a measurement data file:
//// spline spl ("smoothedmeasurements.txt");
//// Temperature field:
//field T("h1");
//// Young’s modulus as a cubic (natural) spline interpolation of the experimental data:
//expression E(spl, T);
//This creates a continuous expression based on discrete data samples. As an application example, if
//measurements of a material stiffness (Young’s modulus E) have been performed for a set of temperatures T (as illustrated in the example above) then this constructor allows to define expression E
//that provides a 3rd order (natural) spline interpolation of Young’s modulus in the measured discrete
//temperature range. Field T can contain any space-dependent temperature profile as long as it is in
//the temperature data range provided.
// Stainless steel 416 BH Curve
// B  - - - H
// 0.000    0.0
// 0.200    318.3
// 0.426    477.5
// 0.761    795.8
// 1.097    1591.6
// 1.233    2387.3
// 1.335    3978.9
// 1.460    7957.8
// 1.590    15915.5
// 1.690    31831.0
// 1.724    44456.3
// 1.740    55704.3

// BH Curve first try
std :: vector<double> Bcurve = {0.0, 0.2,   0.426, 0.761,  1.097, 1.233,  1.335,   1.460,   1.590,   1.690,   1.724,   1.740 , 2.0, 2.5};
std :: vector<double> Hcurve = {0.0, 318.3, 477.5, 795.8, 1591.6, 2387.3, 3978.9, 7957.8, 15915.5, 31831.0, 44456.3, 55704.3 , 58000.0, 60000.0};
spline spl(Bcurve, Hcurve);
spl.write("BHCurve.txt", 5);
cout << "BH Curve..." << endl;
expression BHCurve(spl.getderivative(), norm(curl(az)));
cout << "BH Curve for iron..." << endl;
expression mIron = BHCurve;

    parameter mu;
    mu|all = mu0;
    // Overwrite on non-magnetic regions:
    mu|magnets = mu0;
    mu|iron = mIron;//2000*mu0;
cout << "Formulation" << endl;
    formulation magnetostatics;
    // The strong form of the magnetostatic formulation is curl( 1/mu * curl(a) ) = j, with b = curl(a):
    magnetostatics += integral(all, 1/mu* curl(dof(az)) * curl(tf(az)) );

    // Add the remanent magnetization of the rotor magnet:
    magnetostatics += integral(magnetsdown, -1/mu* minusbremanent * curl(tf(az)));
    magnetostatics += integral(magnetsup, -1/mu* bremanent * curl(tf(az)));
    // magnetostatics += continuitycondition(IronTop, SecondMagnetBottom, az, az, 1, false);

    if(timeit) clk.print("Setup experiment:\t");
    if(autoload == true) {
        cout << "Autoloading cache..." << endl;
        az.loadraw("cache/aaz"+twodigits(alpha)+twodigits(ironAlpha)+".slz.gz", true);
        cout << "Cache loaded" << endl;
        if(timeit) clk.print("CacheLoaded:\t");
    } else {
        cout << "Solving..." << endl;
        magnetostatics.solve();
        az.writeraw(all, "cache/aaz"+twodigits(alpha)+twodigits(ironAlpha)+".slz.gz", true);
        cout << "Solved!" << endl;
        if(timeit) clk.print("Solving time:\t");
    }

    // phiProjection.setdata(all, )


    // expression polito = polarscalar(curl(az));
    // double integration = polito.integrate(firstmagnetup, 4);
    // double integrationdown = polito.integrate(firstmagnetdown, 4);

    // polito.write(bottommagnets, "results/bottommagnets"+twodigits(alpha)+".vtu");
    // cout << "First row INTEGRATION  UP : " << polito.integrate(firstmagnetup, 4) << endl;
    // cout << "First row INTEGRATION DOWN: " << polito.integrate(firstmagnetdown, 4) << endl;
    // cout << " ---------------------------------------------------------------- " << integration << endl;
    // cout << " INTEGRATION  IRON : " << polito.integrate(ironbars, 4) << endl;
    // // cout << " INTEGRATION : " << integrationdown << endl;
    // cout << " ---------------------------------------------------------------- " << integration << endl;
    // cout << " Second row INTEGRATION  UP : " << polito.integrate(secondmagnetup, 4) << endl;
    // cout << " Second row INTEGRATION DOWN: " << polito.integrate(secondmagnetdown, 4) << endl;


    // The magnetostatic force acting on the motor is computed below.
    field magforce("h1xyz");
    magforce.setorder(all, 2);

    // The magnetic force is projected on field 'magforce' on the solid stator region.
    // This is done with a formulation of the type dof*tf - force calculation = 0.
    formulation forceprojection;

    forceprojection += integral(all, dof(magforce)*tf(magforce));
    forceprojection += integral(all, -predefinedmagnetostaticforce(tf(magforce, all), 1/mu*curl(az), mu));

    forceprojection.solve();

    if(timeit) clk.print("Forces calc time:\t");

    expression magforcescalar = polarscalar((expression)magforce);
    // magforcescalar.write(bottomPart, "results/bottomMagForce"+twodigits(alpha)+".vtu");


    // Calculate the torque:
    expression leverarm = array3x1(x,y,0);
    expression torgue = compz(crossproduct(leverarm, magforce));

    double T = torgue.integrate(bottomPart, 4);
    double TBottom = T;
    double TTop = torgue.integrate(topPart, 4);
    double TIron = torgue.integrate(middlesection, 4);
    double TBM = torgue.integrate(bottommagnets, 4);
    double TTM = torgue.integrate(topmagnets, 4);

    double scaleFactor = 1.;
    cout << "--------------------------------------------------------" << endl;
    cout << "TorgueBottom: " << TBottom   * scaleFactor << endl;
    cout << "TorgueIron: " << TIron   * scaleFactor << endl;
    cout << "TorgueTop: " << TTop   * scaleFactor << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------" << endl;

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
    pack->H = curl(az)/mu;
    pack->alpha = alpha;
    pack->ironAlpha = ironAlpha;
    pack->allRegion = selectall();

    // pack->saveData();

    if(timeit) clk.print("Pack Ready:\t");
    return pack;
}

string getFilename(int alpha,int ironAlpha,int argc,char* argv[], bool autoname = false) {
    // FILE MESH 2.2 .msh
    string num = twodigits(alpha);
    string numiron = twodigits(ironAlpha);
    // string test = "models/New.msh";
    string filename;
    if ( argc < 4 && autoname == false ) {
        // take file from args
        filename = (string)"models/" + argv[1];
        cout << "Loading model: " << filename << endl;
    } else {
        filename = (string)"models/parts/" + argv[1] + num + "A" + numiron + ".msh";
        cout << "Loading model(autoname): " << filename << endl;
    }
    return filename;
}

void simulate(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Provide the name of the mesh file!" << endl;
        return;
    }

    string filename = getFilename(0, 0, argc, argv, true);
    mesh mymesh(filename);

    DataPack* pack;
    pack = sparselizard(mymesh, 0, 0, false);
    pack->saveData();
}

int main(int argc, char *argv[])
{
    PetscInitialize(0,{},0,0);
    wallclock clk;

    if(argc < 2) {
        cout << "Need to specify a model name. e.g. ./mcvt Parrot" << endl;
        return 0;
    }

    int start = 0;
    int end = 45;
    int alphaStep = atoi(getVar(argc, argv, 2));
    int pieces = (end-start)/alphaStep;


    string folderName = "July";

    createPath("cache");

    int ironstart = 0;
    int ironend = 45;
    int ironStep = atoi(getVar(argc, argv, 3));
    int ironpieces = (ironend - ironstart)/ironStep;

    DataPack* pack;
    // make only one loop
    for (double ironAlpha = ironstart; ironAlpha <= ironstart+(ironpieces*ironStep); ironAlpha += ironStep)
    {
        DataArray* torguesB = new DataArray(folderName);
        DataArray* torguesI = new DataArray(folderName);
        DataArray* torguesT = new DataArray(folderName);
        for (double alpha = start; alpha <= start+(pieces*alphaStep); alpha += alphaStep)
        {
            cout << "Iron Angle\t" << ironAlpha << "\t Rotor Angle:\t" << alpha << endl;
            string filename = getFilename(alpha, ironAlpha, argc, argv, true);
            mesh mymesh(filename);
            pack = sparselizard(mymesh, alpha, ironAlpha, false, folderName, true);
//            clk.print("Calculations:");
//            cout << "--------------------------------------" << endl << endl;
            pack->saveData();
            torguesB->add(pack->TBottom);
            torguesI->add(pack->TIron);
            torguesT->add(pack->TTop);
            delete pack;
        }
        torguesB->save("TBottom" + twodigits(ironAlpha));
        torguesI->save("TIron" + twodigits(ironAlpha));
        torguesT->save("TTop" + twodigits(ironAlpha));
        delete torguesB, torguesI, torguesT;
    }

    clk.print("Total run time:");
    PetscFinalize();
    return 0;
}

int oldmain(int argc, char *argv[])
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

    string folderName = "December";

    if(argc > 3) {
        alphaStep = atoi(argv[2]);
        ironStep = atoi(argv[3]);
        pieces = atoi(argv[4]);
        ironpieces = atoi(argv[5]);
        folderName = argv[6];
    } else if(argc == 3) {
        folderName = argv[2];
    } else if(argc == 1) {
        cout
        << "Use this software as follows:" << endl
        << "Example command to load only one mesh file" << endl << endl
        << "./loadem.sh FilenameToLoad" << endl << endl
        << "Example command to load only one mesh file save to folderName (default:December)" << endl << endl
        << "./loadem.sh FilenameToLoad FolderName" << endl << endl
        << "Example command to load only multiple mesh file of the name convention: FilenameToLoad{alpha}A{ironAlpha}.msh" << endl << endl << endl
        << "./loadem.sh FilenameToLoad alphaStep ironStep pieces ironpieces" << endl << endl << endl;
        return 0;
    }
    cout << "Settings" << endl
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
            string filename = getFilename(alpha, ironAlpha, argc, argv);
            mesh mymesh(filename);
            // mymesh.scale(regionall(), 0.01,0.01,0.01);
            pack = sparselizard(mymesh, alpha, ironAlpha, false, folderName);
            pack->saveData();
            
            // torgues.push_back(pack->T);
            // torguesBM.push_back(pack->TBottom);
            // torguesIron.push_back(pack->TIron);
            // torguesTM.push_back(pack->TTop);

            // cout << "Iron Angle\t| Mechanical angle [degrees] and torque [Nm]:" << endl;
            // cout << ironAlpha << "\t\t|" << alpha << "\t\t\t\t\t|" << pack->T << endl;
            
            delete pack;
        }

        // writeVec(torguesBM, "results/IronB"+ to_string((int)(ironAlpha))+"A"+ to_string((int)(alphaStep))+"torguesTBottom.csv");
        // writeVec(torguesTM, "results/IronB"+ to_string((int)(ironAlpha))+"A"+ to_string((int)(alphaStep))+"torguesTTop.csv");
        // writeVec(torguesIron, "results/IronB"+ to_string((int)(ironAlpha))+"A"+ to_string((int)(alphaStep))+"torguesTIron.csv");
        
        // torgues.clear();
        // torguesBM.clear();
        // torguesTM.clear();
        // torguesIron.clear();
    }

    clk.print("Total run time:");

    PetscFinalize();
    
    return 0;
}



expression normalize(expression e)
{
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
